#define _CRT_SECURE_NO_WARNINGS
#include "xoshiro256plusplus.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <omp.h>

// Constants
#define STATE_SIZE 16   // 128位
#define BYTE_SIZE 8
#define MAX_LINE_LENGTH 1024
#define MAX_PAIRS (1ULL<<32)
//#define MAX_PAIRS 200000
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
#define round 100
#define AES_TABLE_ROW(OUT_IDX, T0, T1, T2, T3, IN_IDX_0, IN_IDX_1, IN_IDX_2, IN_IDX_3) \
    ciphertext[OUT_IDX] = aes_table[T0][c_temp[IN_IDX_0]] ^ \
                          aes_table[T1][c_temp[IN_IDX_1]] ^ \
                          aes_table[T2][c_temp[IN_IDX_2]] ^ \
                          aes_table[T3][c_temp[IN_IDX_3]]

// Type definitions
typedef uint8_t State[STATE_SIZE];
typedef uint8_t KeySchedule[round + 2][16];

// Structure definitions
typedef struct {
    State first;
    State second;
} StatePair;

// Function declarations
// 函数原型声明部分
void* check_malloc(size_t size);
void print_state(State diff);
void build_lookup_tables(const uint8_t matrix_sol[8], uint8_t matrix_sol_table[256], uint8_t matrix_sbox_table[256], uint8_t xtime_table[256]);
void compute_diff_output(const uint8_t diff_word[16], const uint8_t diff_sol[16], const uint8_t matrix_sol_table[256], State diff_output);
void build_table_aes(const uint8_t matrix_sbox_table[256], const uint8_t time_table[256], uint8_t aes_table[16][256]);
void generate_rcon(uint8_t rcon[round + 2]);
void key_expansion(const uint8_t key[16], const uint8_t matrix_sbox_table[256], const uint8_t rcon[round + 2], KeySchedule round_keys);
static uint64_t splitmix64(uint64_t* x);
void xoshiro256plusplus_init(xoshiro256plusplus_state* state, uint64_t seed, uint64_t thread_id);
static inline uint64_t rotl(const uint64_t x, int k);
uint64_t xoshiro256plusplus_next(xoshiro256plusplus_state* state);
void generate_pair(StatePair* pair, const uint8_t diff_bites[16], xoshiro256plusplus_state* prng);
void encrypt(const uint8_t plaintext[16], State ciphertext, const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const int index[16]);
bool pass_filter(const uint8_t c1[16], const uint8_t c2[16], const uint8_t diff_output[16]);
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint8_t diff_bites[16], const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const uint8_t diff_output[16], const int index[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed);
void guess_subkeys(const uint8_t matrix_sbox_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index);


// Function declarations
void* check_malloc(size_t size) {
    void* ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    return ptr;
}

void print_state(State diff) {
    for (int i = 0; i < 16; ++i) {
        printf("%02X ", diff[i]);
    }
    printf("\n");
}

// 拆出 lookup tables 构建
void build_lookup_tables(const uint8_t matrix_sol[8], uint8_t matrix_sol_table[256], uint8_t matrix_sbox_table[256], uint8_t xtime_table[256]) {
    int i, j;
    for (i = 0; i < 256; i++) {
        uint8_t temp1 = i, temp_value = 0;
        if (i == 0) temp1 = 1;
        else if (i == 1) temp1 = 0;
		
        for (j = 0; j < 8; j++) {
            if (__builtin_parity(matrix_sol[j] & i)) {
                temp_value |= (1 << j);
            }
        }
        matrix_sol_table[i] = temp_value;

        temp_value = 0;
        for (j = 0; j < 8; j++) {
            if (__builtin_parity(matrix_sol[j] & temp1)) {
                temp_value |= (1 << j);
            }
        }
        matrix_sbox_table[i] = temp_value;
    }

    for (i = 0; i < 256; i++) {
        xtime_table[i] = (i << 1) ^ ((i & 0x80) ? 0x1B : 0x00);
    }
}

// 拆出 diff_output 计算函数
void compute_diff_output(const uint8_t diff_word[16], const uint8_t diff_sol[16], const uint8_t matrix_sol_table[256], State diff_output) {
    for (int i = 0; i < 16; ++i) {
        uint8_t diff_temp = diff_sol[i] ^ diff_word[i];
		diff_output[i] = matrix_sol_table[diff_temp];
    }
}

// 表格构建函数
void build_table_aes(const uint8_t matrix_sbox_table[256], const uint8_t time_table[256], uint8_t aes_table[16][256]) {
    for (int i = 0; i < 256; i++) {
        uint8_t x = matrix_sbox_table[i];
        uint8_t tx = time_table[x];
        uint8_t tx_x = tx ^ x;

        // 列混淆
        aes_table[0][i] = tx;
        aes_table[1][i] = x;
        aes_table[2][i] = x;
        aes_table[3][i] = tx_x;
        aes_table[4][i] = tx_x;
        aes_table[5][i] = tx;
        aes_table[6][i] = x;
        aes_table[7][i] = x;
        aes_table[8][i] = x;
        aes_table[9][i] = tx_x;
        aes_table[10][i] = tx;
        aes_table[11][i] = x;
        aes_table[12][i] = x;
        aes_table[13][i] = x;
        aes_table[14][i] = tx_x;
        aes_table[15][i] = tx;
    }
}

// Rcon生成函数
void generate_rcon(uint8_t rcon[round + 2]) {
    rcon[0] = 0x00;  // 占位
    uint8_t value = 0x01;
    for (int i = 1; i <= round + 1; i++) {
        rcon[i] = value;
        value <<= 1;
        if (value & 0x100) {
            value ^= 0x11B;
        }
    }
}

// 密钥扩展函数
void key_expansion(const uint8_t key[16], const uint8_t matrix_sbox_table[256], const uint8_t rcon[round + 2], KeySchedule round_keys) {
    int j;

    // 拷贝初始密钥
    memcpy(round_keys[0], key, 16);

    for (int r = 1; r <= round + 1; r++) {
        uint8_t temp[4];

        // Get last 4 bytes of previous round key
        memcpy(temp, &round_keys[r - 1][12], 4);

        // RotWord
        uint8_t first = temp[0];
        temp[0] = temp[1];
        temp[1] = temp[2];
        temp[2] = temp[3];
        temp[3] = first;

        // SubWord
        temp[0] = matrix_sbox_table[temp[0]];
        temp[1] = matrix_sbox_table[temp[1]];
        temp[2] = matrix_sbox_table[temp[2]];
        temp[3] = matrix_sbox_table[temp[3]];

        // Add Rcon
        temp[0] ^= rcon[r];

        // Generate first 4 bytes
        round_keys[r][0] = round_keys[r - 1][0] ^ temp[0];
        round_keys[r][1] = round_keys[r - 1][1] ^ temp[1];
        round_keys[r][2] = round_keys[r - 1][2] ^ temp[2];
        round_keys[r][3] = round_keys[r - 1][3] ^ temp[3];

        // Generate remaining 12 bytes
        for (j = 4; j < 16; ++j) {
            round_keys[r][j] = round_keys[r - 1][j] ^ round_keys[r][j - 4];
        }
    }
}

// 加速splitmix64
static uint64_t splitmix64(uint64_t* x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// xoshiro256plusplus密钥编排算法
void xoshiro256plusplus_init(xoshiro256plusplus_state* state, uint64_t seed, uint64_t thread_id) {
    uint64_t sm64 = seed ^ (thread_id * 0xA3C59AC2ULL);
    for (int i = 0; i < 4; ++i)
        state->s[i] = splitmix64(&sm64);
}

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

// xoshiro256plusplus随机数生成函数
uint64_t xoshiro256plusplus_next(xoshiro256plusplus_state* state) {
    uint64_t* s = state->s;
    uint64_t result = rotl(s[0] + s[3], 23) + s[0];

    uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
}

// 生成随机明文结构
void generate_pair(StatePair* pair, const uint8_t diff_bites[16], xoshiro256plusplus_state* prng) {  
    uint64_t r1 = xoshiro256plusplus_next(prng);
    uint64_t r2 = xoshiro256plusplus_next(prng);
    // 拆分成16字节
    pair->first[0] = (r1 >> 0) & 0xFF;
    pair->first[1] = (r1 >> 8) & 0xFF;
    pair->first[2] = (r1 >> 16) & 0xFF;
    pair->first[3] = (r1 >> 24) & 0xFF;
    pair->first[4] = (r1 >> 32) & 0xFF;
    pair->first[5] = (r1 >> 40) & 0xFF;
    pair->first[6] = (r1 >> 48) & 0xFF;
    pair->first[7] = (r1 >> 56) & 0xFF;
    pair->first[8] = (r2 >> 0) & 0xFF;
    pair->first[9] = (r2 >> 8) & 0xFF;
    pair->first[10] = (r2 >> 16) & 0xFF;
    pair->first[11] = (r2 >> 24) & 0xFF;
    pair->first[12] = (r2 >> 32) & 0xFF;
    pair->first[13] = (r2 >> 40) & 0xFF;
    pair->first[14] = (r2 >> 48) & 0xFF;
    pair->first[15] = (r2 >> 56) & 0xFF;
    
    for (int k = 0; k < 16; ++k) {
        pair->second[k] = pair->first[k] ^ diff_bites[k];
    }
}

// round轮aes加密
void encrypt(const uint8_t plaintext[16], State ciphertext, const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const int index[16]) {
    memcpy(ciphertext, plaintext, STATE_SIZE);
    // XOR with initial key
    int i, j;

    for (j = 0; j < 16; ++j) {
        ciphertext[j] ^= key_schedule[0][j];
    }

    // Perform AES rounds
    for (int r = 0; r < round; ++r) {
        State c_temp;
        memcpy(c_temp, ciphertext, STATE_SIZE);

        AES_TABLE_ROW(0, 0, 4, 8, 12, 0, 5, 10, 15);
        AES_TABLE_ROW(1, 1, 5, 9, 13, 0, 5, 10, 15);
        AES_TABLE_ROW(2, 2, 6, 10, 14, 0, 5, 10, 15);
        AES_TABLE_ROW(3, 3, 7, 11, 15, 0, 5, 10, 15);
        AES_TABLE_ROW(4, 0, 4, 8, 12, 4, 9, 14, 3);
        AES_TABLE_ROW(5, 1, 5, 9, 13, 4, 9, 14, 3);
        AES_TABLE_ROW(6, 2, 6, 10, 14, 4, 9, 14, 3);
        AES_TABLE_ROW(7, 3, 7, 11, 15, 4, 9, 14, 3);
        AES_TABLE_ROW(8, 0, 4, 8, 12, 8, 13, 2, 7);
        AES_TABLE_ROW(9, 1, 5, 9, 13, 8, 13, 2, 7);
        AES_TABLE_ROW(10, 2, 6, 10, 14, 8, 13, 2, 7);
        AES_TABLE_ROW(11, 3, 7, 11, 15, 8, 13, 2, 7);
        AES_TABLE_ROW(12, 0, 4, 8, 12, 12, 1, 6, 11);
        AES_TABLE_ROW(13, 1, 5, 9, 13, 12, 1, 6, 11);
        AES_TABLE_ROW(14, 2, 6, 10, 14, 12, 1, 6, 11);
        AES_TABLE_ROW(15, 3, 7, 11, 15, 12, 1, 6, 11);

        for (j = 0; j < 16; ++j) {
            ciphertext[j] ^= key_schedule[r + 1][j];
        }
    }

    // Final round
	uint8_t temp[16];
    for (i = 0; i < 16; ++i) {
        temp[i] = ciphertext[index[i]];
    }
    for (i = 0; i < 16; ++i) {
        ciphertext[i] = matrix_sbox_table[temp[i]];
	}
    for (j = 0; j < 16; ++j) {
        ciphertext[j] ^= key_schedule[round + 1][j];
    }
}

// 判断一对加密后是否满足差分
bool pass_filter(const uint8_t c1[16], const uint8_t c2[16], const uint8_t diff_output[16]) {
    for (int i = 0; i < STATE_SIZE; ++i) {
        if ((c1[i] ^ c2[i]) != diff_output[i]) return false;
    }
    return true;
}

// 单独将加密过程和筛选密文过程分离出来
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint8_t diff_bites[16], const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const uint8_t diff_output[16], const int index[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed) {
    size_t i, j;
    // 限制线程数为物理核心数
    int max_threads = omp_get_max_threads();
    if (max_threads > 32) max_threads = 32;  // 调整为你的服务器物理核心数

    // 使用更大的线程本地buffer
    size_t thread_bufsize = (max_filtered_count / max_threads) + 1024;
    StatePair** thread_buffers = malloc(max_threads * sizeof(StatePair*));
    size_t* thread_counts = calloc(max_threads, sizeof(size_t));

    for (int t = 0; t < max_threads; ++t) {
        thread_buffers[t] = malloc(thread_bufsize * sizeof(StatePair));
    }

    // 推荐大BATCH_SIZE，现代主机16~128k都可以
    const size_t BATCH_SIZE = 1024; // 可根据内存调优

    // 用于每个线程记录其负责的全局pair区间
    size_t pairs_per_thread = (total_pairs + max_threads - 1) / max_threads;

    #pragma omp parallel num_threads(max_threads)
    {
        int tid = omp_get_thread_num();
        size_t my_start = tid * pairs_per_thread;
        size_t my_end = my_start + pairs_per_thread;
        if (my_end > total_pairs) my_end = total_pairs;

        StatePair* local_buf = thread_buffers[tid];
        size_t local_count = 0;

        // 初始化线程本地PRNG
        xoshiro256plusplus_state prng;
        xoshiro256plusplus_init(&prng, global_seed, tid);

        // 批处理临时数组
        StatePair batch[BATCH_SIZE];
        State c1[BATCH_SIZE], c2[BATCH_SIZE];

        //#pragma omp for schedule(dynamic, 1024)
        for (i = my_start; i < my_end; i += BATCH_SIZE) {
            size_t batch_count = (i + BATCH_SIZE > my_end) ? (my_end - i) : BATCH_SIZE;

            // 1. 批量生成明文对
            for (j = 0; j < batch_count; ++j) {
                // 每次用prng生成128位明文
                uint64_t r1 = xoshiro256plusplus_next(&prng);
                uint64_t r2 = xoshiro256plusplus_next(&prng);
                for (int k = 0; k < 8; ++k) {
                    batch[j].first[k] = (r1 >> (8 * k)) & 0xFF;
                    batch[j].first[k + 8] = (r2 >> (8 * k)) & 0xFF;
                    batch[j].second[k] = batch[j].first[k] ^ diff_bites[k];
                    batch[j].second[k + 8] = batch[j].first[k + 8] ^ diff_bites[k + 8];
                }
            }

            // 2. 批量加密
            for (j = 0; j < batch_count; ++j) {
                encrypt(batch[j].first, c1[j], key_schedule, aes_table, matrix_sbox_table, index);
                encrypt(batch[j].second, c2[j], key_schedule, aes_table, matrix_sbox_table, index);
            }

            // 3. 批量过滤
            for (j = 0; j < batch_count; ++j) {
                if (pass_filter(c1[j], c2[j], diff_output)) {
                    if (local_count < thread_bufsize) {
                        memcpy(local_buf[local_count].first, c1[j], STATE_SIZE);
                        memcpy(local_buf[local_count].second, c2[j], STATE_SIZE);
                        local_count++;
                    }
                }
            }
        }
        thread_counts[tid] = local_count;
    }

    // 合并结果
    *filtered_count = 0;
    for (int t = 0; t < max_threads; ++t) {
        size_t copy_count = thread_counts[t];
        if (*filtered_count + copy_count > max_filtered_count) {
            copy_count = max_filtered_count - *filtered_count;
        }
        if (copy_count > 0) {
            memcpy(&filtered_pairs[*filtered_count],
                thread_buffers[t],
                copy_count * sizeof(StatePair));
            *filtered_count += copy_count;
        }
        free(thread_buffers[t]);
    }
    free(thread_buffers);
    free(thread_counts);
}

// 单独将猜测子密钥过程分离出来
void guess_subkeys(const uint8_t matrix_sbox_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index) {
    int i, j;
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    uint8_t diff_det1 = diff_output[index_inv[index_diff_word[0]]] ^ 1;
    uint8_t diff_det2 = diff_output[index_inv[index_diff_word[1]]] ^ 1;

    int* cand_key_match = (int*)calloc(1 << 16, sizeof(int));

    for (i = 0; i < (1 << 16); ++i) {
        uint8_t key1 = i & 0xFF;
        uint8_t key2 = (i >> 8) & 0xFF;
        int num = 0;

        for (j = 0; j < pairs_cand_count; ++j) {
            uint8_t ct1, ct2, ct3, ct4;

            ct1 = pairs_cand[j].first[index_inv[index_diff_word[0]]] ^ key1;
            ct2 = pairs_cand[j].first[index_inv[index_diff_word[1]]] ^ key2;
            ct3 = pairs_cand[j].second[index_inv[index_diff_word[0]]] ^ key1;
            ct4 = pairs_cand[j].second[index_inv[index_diff_word[1]]] ^ key2;

            
            ct1 = matrix_sbox_table[ct1];
            ct2 = matrix_sbox_table[ct2];
            ct3 = matrix_sbox_table[ct3];
            ct4 = matrix_sbox_table[ct4];

            bool match = true;
            if ((ct1 ^ ct3) != diff_det1) match = false;
            if ((ct2 ^ ct4) != diff_det2) match = false;

            if (match) num++;
        }
        cand_key_match[i] = num;
    }

    for (i = 1; i < (1 << 16); ++i) {
        if (cand_key_match[i] > *max_val) {
            *max_val = cand_key_match[i];
        }
    }
    for (i = 0; i < (1 << 16); ++i) {
        if (cand_key_match[i] == *max_val) {
            max_index[*max_index_count] = i;
            (*max_index_count)++;
        }
    }

    printf("Maximum element: %d\n", *max_val);
    printf("Maximum index count: %d\n", *max_index_count);
    for (i = 0; i < *max_index_count; ++i) {
        printf("Index of maximum element: %d\n", max_index[i]);

        uint8_t cand_key1 = max_index[i] & 0xFF;
        uint8_t cand_key2 = (max_index[i] >> 8) & 0xFF;
        printf("Candidate key1: 0x%02X, Candidate key2: 0x%02X\n", cand_key1, cand_key2);
    }

    free(cand_key_match);
}

int main() {
    // 2. 初始化变量和数据结构
    int i, j;
    srand((unsigned int)time(NULL));
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    State diff_sol = { 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00, 0x00 };
    State diff_word = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    uint8_t matrix_sol[] = { 0x97, 0x92, 0x94, 0x98, 0x10, 0x20, 0x40, 0x80 };
    uint64_t global_seed = 0xC7E6A3B2F19D4E58ULL; // 全局种子，用于xoshiro256plusplus初始化

    // 4. 构建查找表
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t matrix_sbox_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, matrix_sbox_table, xtime_table);

    // 3. 计算 diff_output
    State diff_output;
    compute_diff_output(diff_word, diff_sol, matrix_sol_table, diff_output);
    print_state(diff_output);

    // 5. 构建AES表
    uint8_t aes_table[16][256] = { { 0x00 } };

    build_table_aes(matrix_sbox_table, xtime_table, aes_table);

    // 6. 生成主密钥以及密钥拓展
    // 生成轮常数rcon
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sbox_table, rcon, key_schedule);

    // 7. 生成明文 & 加密过程 & 筛选 pairs
    StatePair* filtered_pairs = malloc(MAX_FILTERED_PAIRS * sizeof(StatePair));
    size_t filtered_count = 0;

    clock_t start = clock();

    generate_encrypt_and_filter_stream(MAX_PAIRS, diff_sol, key_schedule, aes_table, matrix_sbox_table, diff_output, index, filtered_pairs, &filtered_count, MAX_FILTERED_PAIRS, global_seed);
    printf("最终通过筛选的明密文对数：%zu\n", filtered_count);
    clock_t end = clock();
    double duration = ((double)(end - start)) / CLOCKS_PER_SEC * 1000.0;
    printf("循环耗时: %.2f 毫秒\n", duration);

    // 8. 猜测密钥
    int index_diff_word[2];
    int index_diff_word_count = 0;
    for (int i = 0; i < 16; ++i) {
        if (diff_word[i]) {
            index_diff_word[index_diff_word_count++] = i;
        }
    }
    int max_val = 0;
    int max_index[1 << 16];
    int max_index_count = 0;

	guess_subkeys(matrix_sbox_table, filtered_pairs, filtered_count, diff_output, index_diff_word, &max_val, &max_index_count, max_index);

    //uint8_t right_key1 = key_schedule[round + 1][index_inv[index_diff_word[0]]];
    //uint8_t right_key2 = key_schedule[round + 1][index_inv[index_diff_word[1]]];
    //printf("Right key1: 0x%02X, Right key2: 0x%02X\n", right_key1, right_key2); // 大写十六进制输出：0xAB
    /*for (int i = 0; i < round + 2; ++i) {
        for (int j = 0; j < 16; ++j) {
            printf("0x%02X, ", key_schedule[i][j]);
        }
        printf("\n");
    }
    printf("\n");*/

    // 10. 导出全部候选子密钥
    FILE* fp = fopen("cand_key_0.txt", "w");   // 打开文件，"w" 表示写入模式

    // 将数组元素逐个写入文件
    for (i = 0; i < max_index_count; i++) {
        fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
    }

    fclose(fp);  // 关闭文件
    printf("数组已成功写入到 cand_key_0.txt\n");

    // 11. 释放内存
    free(filtered_pairs);

    return 0;
}