#define _CRT_SECURE_NO_WARNINGS
#include "xoshiro256plusplus.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <mpi.h>

// Constants
#define STATE_SIZE 16   // 128位
#define BYTE_SIZE 8
#define BATCH_SIZE 1024  // 将批处理大小定义为常量
#define MAX_LINE_LENGTH 1024
//#define MAX_PAIRS (1ULL<<45)
#define MAX_PAIRS 500000
#define round 21
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
#define AES_TTABLE_ROW(out, t0, t1, t2, t3, rk) \
    (T0[(t0)] ^ T1[(t1)] ^ T2[(t2)] ^ T3[(t3)] ^ (rk))

// Type definitions
typedef uint8_t State[16] __attribute__((aligned(8)));
typedef uint8_t KeySchedule[round + 2][16] __attribute__((aligned(8)));

// Structure definitions
typedef struct __attribute__((aligned(16))) {
    uint8_t first[16];
    uint8_t second[16];
} StatePair;

// 错误处理宏
#define MPI_CHECK(call) \
    do { \
        int err = call; \
        if (err != MPI_SUCCESS) { \
            char error_string[MPI_MAX_ERROR_STRING]; \
            int length; \
            MPI_Error_string(err, error_string, &length); \
            fprintf(stderr, "MPI error at %s:%d: %s\n", __FILE__, __LINE__, error_string); \
            MPI_Abort(MPI_COMM_WORLD, err); \
        } \
    } while (0)

// Function declarations
// 函数原型声明部分
void print_state(State diff);
void build_lookup_tables(const uint8_t matrix_sol[8], uint8_t matrix_sol_table[256], uint8_t matrix_sbox_table[256], uint8_t xtime_table[256]);
void build_table_aes(const uint8_t matrix_sbox_table[256], const uint8_t time_table[256]);
void generate_rcon(uint8_t rcon[round + 2]);
void key_expansion(const uint8_t key[16], const uint8_t matrix_sbox_table[256], const uint8_t rcon[round + 2], KeySchedule round_keys);
static inline uint64_t splitmix64(uint64_t* x);
void xoshiro256plusplus_init(xoshiro256plusplus_state* state, uint64_t seed, uint64_t thread_id);
static inline uint64_t rotl(const uint64_t x, int k);
uint64_t xoshiro256plusplus_next(xoshiro256plusplus_state* state);
void encrypt(const uint8_t plaintext[16], uint8_t ciphertext[16], const uint32_t key_schedule_32[round + 2][4], const KeySchedule key_schedule, const uint8_t matrix_sbox_table[256]);
bool pass_filter(const uint8_t c1[16], const uint8_t c2[16], const uint64_t diff_output_64[2]);
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint64_t diff_bites[2], const uint32_t key_schedule_32[round + 2][4], const KeySchedule key_schedule, const uint8_t matrix_sbox_table[256], const uint64_t diff_output_64[2], StatePair* batch, StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed, int mpi_rank, int mpi_size);
void guess_subkeys(const uint8_t matrix_sbox_table[256], const StatePair* pairs_cand, int pairs_cand_count, const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index);

// 全局变量定义（main函数之外）
uint32_t T0[256], T1[256], T2[256], T3[256];
int row_idx[16] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
int row_idx_inv[16] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
uint8_t diff_output[16] = { 0x00, 0x00, 0x4e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4e, 0x00, 0x00, 0x00, 0x00, 0x00 };

void print_state(State diff) {
    for (int i = 0; i < 16; ++i) {
        printf("%02X ", diff[i]);
    }
    printf("\n");
}

// 拆出 lookup tables 构建
void build_lookup_tables(const uint8_t matrix_sol[8], uint8_t matrix_sol_table[256], uint8_t matrix_sbox_table[256], uint8_t xtime_table[256]) {
    int i, j;
    // 1. 先预算基向量的映射
    uint8_t base[8];
	uint8_t temp1, temp_value;
    for (i = 0; i < 8; ++i) {
        // 只设置第k位，其它为0
        temp1 = (1 << i);
        temp_value = 0;
        for (j = 0; j < 8; ++j) {
            if (__builtin_parity(matrix_sol[j] & temp1)) {
                temp_value |= (1 << j);
            }
        }
        base[i] = temp_value; // base[i] = matrix_sol * (1 << i)
    }

    // 2. 用基向量异或展开所有x
    for (i = 0; i < 256; ++i) {
        temp_value = 0;
        for (j = 0; j < 8; ++j) {
            if (i & (1 << j)) {
                temp_value ^= base[j];
            }
        }
        matrix_sol_table[i] = temp_value;
    }

    // 3. matrix_sbox_table[x] = matrix_sol_table[simple_swap(x)]
    // 你原来是 i==0 <-> 1，其它不变
    matrix_sbox_table[0] = matrix_sol_table[1];
    matrix_sbox_table[1] = matrix_sol_table[0];
    for (i = 2; i < 256; ++i) {
        matrix_sbox_table[i] = matrix_sol_table[i];
    }

    // 4. xtime_table
    for (i = 0; i < 256; ++i) {
        xtime_table[i] = (i << 1) ^ ((i & 0x80) ? 0x1B : 0x00);
    }
}

// 表格构建函数
void build_table_aes(const uint8_t matrix_sbox_table[256], const uint8_t time_table[256]) {
    int i;
    uint8_t x, tx, tx_x;
    for (i = 0; i < 256; ++i) {
        x = matrix_sbox_table[i];
        tx = time_table[x];
        tx_x = tx ^ x;
        // 你需要根据你原来的列混淆规则分别合成T0~T3
        // 下面是示例，需根据你的实际列混淆顺序调整
        T0[i] = (tx << 24) | (x << 16) | (x << 8) | tx_x;
        T1[i] = (tx_x << 24) | (tx << 16) | (x << 8) | x;
        T2[i] = (x << 24) | (tx_x << 16) | (tx << 8) | x;
        T3[i] = (x << 24) | (x << 16) | (tx_x << 8) | tx;
    }
}

// Rcon生成函数
void generate_rcon(uint8_t rcon[round + 2]) {
    rcon[0] = 0x00;  // 占位
    uint8_t value = 0x01;
    int i;
    for (i = 1; i <= round + 1; i++) {
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
    int r;

    for (r = 1; r <= round + 1; r++) {
        round_keys[r][0] = matrix_sbox_table[round_keys[r - 1][13]];
        round_keys[r][1] = matrix_sbox_table[round_keys[r - 1][14]];
        round_keys[r][2] = matrix_sbox_table[round_keys[r - 1][15]];
        round_keys[r][3] = matrix_sbox_table[round_keys[r - 1][12]];
        round_keys[r][0] ^= rcon[r];

        // Generate first 4 bytes
        round_keys[r][0] ^= round_keys[r - 1][0];
        round_keys[r][1] ^= round_keys[r - 1][1];
        round_keys[r][2] ^= round_keys[r - 1][2];
        round_keys[r][3] ^= round_keys[r - 1][3];

        // Generate remaining 12 bytes
        round_keys[r][4] = round_keys[r][0] ^ round_keys[r - 1][4];
        round_keys[r][5] = round_keys[r][1] ^ round_keys[r - 1][5];
        round_keys[r][6] = round_keys[r][2] ^ round_keys[r - 1][6];
        round_keys[r][7] = round_keys[r][3] ^ round_keys[r - 1][7];
        round_keys[r][8] = round_keys[r][4] ^ round_keys[r - 1][8];
        round_keys[r][9] = round_keys[r][5] ^ round_keys[r - 1][9];
        round_keys[r][10] = round_keys[r][6] ^ round_keys[r - 1][10];
        round_keys[r][11] = round_keys[r][7] ^ round_keys[r - 1][11];
        round_keys[r][12] = round_keys[r][8] ^ round_keys[r - 1][12];
        round_keys[r][13] = round_keys[r][9] ^ round_keys[r - 1][13];
        round_keys[r][14] = round_keys[r][10] ^ round_keys[r - 1][14];
        round_keys[r][15] = round_keys[r][11] ^ round_keys[r - 1][15];
    }
}

// 加速splitmix64
static inline uint64_t splitmix64(uint64_t* x) {
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

// round轮aes加密
// 假设State为uint8_t[16]，KeySchedule为uint8_t[rounds+2][16]
void encrypt(const uint8_t plaintext[16], uint8_t ciphertext[16], const uint32_t key_schedule_32[round + 2][4], const KeySchedule key_schedule, const uint8_t matrix_sbox_table[256]) {
    int i, r;
    // 打包明文并加初始密钥
    uint32_t* s = (uint32_t*)plaintext; // 直接在明文buffer上操作;
    s[0] = __builtin_bswap32(s[0]);
    s[1] = __builtin_bswap32(s[1]);
    s[2] = __builtin_bswap32(s[2]);
    s[3] = __builtin_bswap32(s[3]);
	s[0] ^= key_schedule_32[0][0];
    s[1] ^= key_schedule_32[0][1];
    s[2] ^= key_schedule_32[0][2];
    s[3] ^= key_schedule_32[0][3];

    uint32_t t[4];
    for (int r = 1; r < round + 1; ++r) {
        t[0] = T0[(s[0] >> 24) & 0xFF] ^ T1[(s[1] >> 16) & 0xFF] ^ T2[(s[2] >> 8) & 0xFF] ^ T3[s[3] & 0xFF] ^ key_schedule_32[r][0];
        t[1] = T0[(s[1] >> 24) & 0xFF] ^ T1[(s[2] >> 16) & 0xFF] ^ T2[(s[3] >> 8) & 0xFF] ^ T3[s[0] & 0xFF] ^ key_schedule_32[r][1];
        t[2] = T0[(s[2] >> 24) & 0xFF] ^ T1[(s[3] >> 16) & 0xFF] ^ T2[(s[0] >> 8) & 0xFF] ^ T3[s[1] & 0xFF] ^ key_schedule_32[r][2];
        t[3] = T0[(s[3] >> 24) & 0xFF] ^ T1[(s[0] >> 16) & 0xFF] ^ T2[(s[1] >> 8) & 0xFF] ^ T3[s[2] & 0xFF] ^ key_schedule_32[r][3];
        // 直接写回明文buffer
        s[0] = t[0];
        s[1] = t[1];
        s[2] = t[2];
        s[3] = t[3];
    }
    // 最后一轮（ShiftRows+S盒+加密钥）
    // 行移位
    uint8_t tmp[16];
    tmp[row_idx_inv[0]] = (s[0] >> 24) & 0xFF;
    tmp[row_idx_inv[1]] = (s[0] >> 16) & 0xFF;
    tmp[row_idx_inv[2]] = (s[0] >> 8) & 0xFF;
    tmp[row_idx_inv[3]] = s[0] & 0xFF;
    tmp[row_idx_inv[4]] = (s[1] >> 24) & 0xFF;
    tmp[row_idx_inv[5]] = (s[1] >> 16) & 0xFF;
    tmp[row_idx_inv[6]] = (s[1] >> 8) & 0xFF;
    tmp[row_idx_inv[7]] = s[1] & 0xFF;
    tmp[row_idx_inv[8]] = (s[2] >> 24) & 0xFF;
    tmp[row_idx_inv[9]] = (s[2] >> 16) & 0xFF;
    tmp[row_idx_inv[10]] = (s[2] >> 8) & 0xFF;
    tmp[row_idx_inv[11]] = s[2] & 0xFF;
    tmp[row_idx_inv[12]] = (s[3] >> 24) & 0xFF;
    tmp[row_idx_inv[13]] = (s[3] >> 16) & 0xFF;
    tmp[row_idx_inv[14]] = (s[3] >> 8) & 0xFF;
    tmp[row_idx_inv[15]] = s[3] & 0xFF;
    
    // S盒+加密钥
    ciphertext[0] = matrix_sbox_table[tmp[0]] ^ key_schedule[round + 1][0];
    ciphertext[1] = matrix_sbox_table[tmp[1]] ^ key_schedule[round + 1][1];
    ciphertext[2] = matrix_sbox_table[tmp[2]] ^ key_schedule[round + 1][2];
    ciphertext[3] = matrix_sbox_table[tmp[3]] ^ key_schedule[round + 1][3];
    ciphertext[4] = matrix_sbox_table[tmp[4]] ^ key_schedule[round + 1][4];
    ciphertext[5] = matrix_sbox_table[tmp[5]] ^ key_schedule[round + 1][5];
    ciphertext[6] = matrix_sbox_table[tmp[6]] ^ key_schedule[round + 1][6];
    ciphertext[7] = matrix_sbox_table[tmp[7]] ^ key_schedule[round + 1][7];
    ciphertext[8] = matrix_sbox_table[tmp[8]] ^ key_schedule[round + 1][8];
    ciphertext[9] = matrix_sbox_table[tmp[9]] ^ key_schedule[round + 1][9];
    ciphertext[10] = matrix_sbox_table[tmp[10]] ^ key_schedule[round + 1][10];
    ciphertext[11] = matrix_sbox_table[tmp[11]] ^ key_schedule[round + 1][11];
    ciphertext[12] = matrix_sbox_table[tmp[12]] ^ key_schedule[round + 1][12];
    ciphertext[13] = matrix_sbox_table[tmp[13]] ^ key_schedule[round + 1][13];
    ciphertext[14] = matrix_sbox_table[tmp[14]] ^ key_schedule[round + 1][14];
    ciphertext[15] = matrix_sbox_table[tmp[15]] ^ key_schedule[round + 1][15];
}

// 判断一对加密后是否满足差分
bool pass_filter(const uint8_t c1[16], const uint8_t c2[16], const uint64_t diff_output_64[2]) {
    const uint64_t* a = (const uint64_t*)c1;
    const uint64_t* b = (const uint64_t*)c2;
    return ((a[0] ^ b[0]) == diff_output_64[0]) && ((a[1] ^ b[1]) == diff_output_64[1]);
}

// 单独将加密过程和筛选密文过程分离出来
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint64_t diff_bites[2], const uint32_t key_schedule_32[round + 2][4], const KeySchedule key_schedule, const uint8_t matrix_sbox_table[256], const uint64_t diff_output_64[2], StatePair* batch, StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed, int mpi_rank, int mpi_size) {
    size_t i, j;
    // 1. 计算本进程需要处理的pair区间
    size_t pairs_per_rank = (total_pairs + mpi_size - 1) / mpi_size;
    size_t my_start = mpi_rank * pairs_per_rank;
    size_t my_end = my_start + pairs_per_rank;
    if (my_end > total_pairs) my_end = total_pairs;

    // 用于记录本地筛选结果的计数器
    size_t local_count = 0;

    // 3. 初始化PRNG
    xoshiro256plusplus_state prng;
    xoshiro256plusplus_init(&prng, global_seed, mpi_rank);

    size_t batch_count;
    uint64_t r1, r2;
    int k;
    // 4. 主循环：批量生成、加密、筛选
    for (i = my_start; i < my_end; i += BATCH_SIZE) {
        batch_count = (i + BATCH_SIZE > my_end) ? (my_end - i) : BATCH_SIZE;

        // 4.1. 批量生成明文对
        for (j = 0; j < batch_count; ++j) {
            // 每次用prng生成128位明文
            *(uint64_t*)(batch[j].first) = xoshiro256plusplus_next(&prng);
            *(uint64_t*)(batch[j].first + 8) = xoshiro256plusplus_next(&prng);

            ((uint64_t*)batch[j].second)[0] = ((uint64_t*)batch[j].first)[0] ^ diff_bites[0];
            ((uint64_t*)batch[j].second)[1] = ((uint64_t*)batch[j].first)[1] ^ diff_bites[1];
        }

        // 4.2. 批量加密和过滤
        for (j = 0; j < batch_count; ++j) {
            encrypt(batch[j].first, batch[j].first, key_schedule_32, key_schedule, matrix_sbox_table);
            encrypt(batch[j].second, batch[j].second, key_schedule_32, key_schedule, matrix_sbox_table);

            if (pass_filter(batch[j].first, batch[j].second, diff_output_64)) {
                if (local_count < max_filtered_count) {  // 使用传入的max_filtered_count
                    // 代替两次 memcpy
                    ((uint64_t*)filtered_pairs[local_count].first)[0] = ((uint64_t*)batch[j].first)[0];
                    ((uint64_t*)filtered_pairs[local_count].first)[1] = ((uint64_t*)batch[j].first)[1];
                    ((uint64_t*)filtered_pairs[local_count].second)[0] = ((uint64_t*)batch[j].second)[0];
                    ((uint64_t*)filtered_pairs[local_count].second)[1] = ((uint64_t*)batch[j].second)[1];

                    local_count++;
                }
            }
        }
    }

    // 5. 返回筛选结果数量
    *filtered_count = local_count;
}

// 单独将猜测子密钥过程分离出来
void guess_subkeys(const uint8_t matrix_sbox_table[256], const StatePair* pairs_cand, int pairs_cand_count, const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index) {
    int i, j, num;
    uint8_t diff_det1 = diff_output[row_idx_inv[index_diff_word[0]]] ^ 1;
    uint8_t diff_det2 = diff_output[row_idx_inv[index_diff_word[1]]] ^ 1;
    uint8_t key1, key2, ct1, ct2, ct3, ct4;
    bool match;

    int* cand_key_match = (int*)calloc(1 << 16, sizeof(int));

    for (i = 0; i < (1 << 16); ++i) {
        key1 = i & 0xFF;
        key2 = (i >> 8) & 0xFF;
        num = 0;

        for (j = 0; j < pairs_cand_count; ++j) {
            ct1, ct2, ct3, ct4;

            ct1 = pairs_cand[j].first[row_idx_inv[index_diff_word[0]]] ^ key1;
            ct2 = pairs_cand[j].first[row_idx_inv[index_diff_word[1]]] ^ key2;
            ct3 = pairs_cand[j].second[row_idx_inv[index_diff_word[0]]] ^ key1;
            ct4 = pairs_cand[j].second[row_idx_inv[index_diff_word[1]]] ^ key2;

            ct1 = matrix_sbox_table[ct1];
            ct2 = matrix_sbox_table[ct2];
            ct3 = matrix_sbox_table[ct3];
            ct4 = matrix_sbox_table[ct4];

            match = true;
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

	uint8_t cand_key1, cand_key2;
    printf("Maximum element: %d\n", *max_val);
    printf("Maximum index count: %d\n", *max_index_count);
    for (i = 0; i < *max_index_count; ++i) {
        printf("Index of maximum element: %d\n", max_index[i]);

        cand_key1 = max_index[i] & 0xFF;
        cand_key2 = (max_index[i] >> 8) & 0xFF;
        printf("Candidate key1: 0x%02X, Candidate key2: 0x%02X\n", cand_key1, cand_key2);
    }

    free(cand_key_match);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Total MPI processes: %d\n\n", size);
    }

    // 1. 初始化变量和数据结构
    int i, j;
    srand((unsigned int)time(NULL));
    State diff_sol = { 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00, 0x00 };
    State diff_plain = { 0x26, 0x00, 0xf7, 0x00, 0x00, 0x9f, 0x00, 0x01, 0x26, 0x00, 0xf7, 0x00, 0x00, 0x9f, 0x00, 0x01 };
    uint64_t diff_plain_64[2];
    memcpy(diff_plain_64, diff_plain, 16);
    State diff_word = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    uint8_t matrix_sol[] = { 0x97, 0x92, 0x94, 0x98, 0x10, 0x20, 0x40, 0x80 };
    uint64_t global_seed = 0xC7E6A3B2F19D4E58ULL; // 全局种子，用于xoshiro256plusplus初始化
    uint64_t diff_output_64[2]; // 用于64位比较
    memcpy(diff_output_64, diff_output, 16);

    // 2. 构建查找表
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t matrix_sbox_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, matrix_sbox_table, xtime_table);

    uint8_t aes_table[16][256] = { { 0x00 } };
    build_table_aes(matrix_sbox_table, xtime_table);
    
    if (rank == 0) {
        print_state(diff_output);
    }

    // 4. 生成主密钥以及密钥拓展
    // 生成轮常数rcon
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

    uint8_t key[16] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sbox_table, rcon, key_schedule);
    uint32_t key_schedule_32[round + 2][4];
    for (i = 0; i < round + 2; ++i) {
		key_schedule_32[i][0] = (key_schedule[i][0] << 24) | (key_schedule[i][1] << 16) | (key_schedule[i][2] << 8) | key_schedule[i][3];
		key_schedule_32[i][1] = (key_schedule[i][4] << 24) | (key_schedule[i][5] << 16) | (key_schedule[i][6] << 8) | key_schedule[i][7];
        key_schedule_32[i][2] = (key_schedule[i][8] << 24) | (key_schedule[i][9] << 16) | (key_schedule[i][10] << 8) | key_schedule[i][11];
        key_schedule_32[i][3] = (key_schedule[i][12] << 24) | (key_schedule[i][13] << 16) | (key_schedule[i][14] << 8) | key_schedule[i][15];
    }

    // 5. 各进程独立生成、加密、筛选
    // 创建MPI派生类型
    MPI_Datatype pair_type;
    MPI_Type_contiguous(32, MPI_BYTE, &pair_type);
    MPI_Type_commit(&pair_type);

    MPI_Barrier(MPI_COMM_WORLD); // 所有进程同步，开始计时
    double t_start = MPI_Wtime();

    StatePair* local_pairs = malloc(MAX_FILTERED_PAIRS * sizeof(StatePair));
    size_t local_count = 0;

    // 2. 分配对齐的本地缓冲区
    StatePair* batch = (StatePair*)aligned_alloc(64, BATCH_SIZE * sizeof(StatePair));

    generate_encrypt_and_filter_stream(MAX_PAIRS, diff_plain_64, key_schedule_32, key_schedule, matrix_sbox_table, diff_output_64, batch, local_pairs, &local_count, MAX_FILTERED_PAIRS, global_seed, rank, size);

    MPI_Barrier(MPI_COMM_WORLD); // 所有进程结束主计算
    double t_end = MPI_Wtime();
    if (rank == 0) {
        printf("主计算部分耗时: %.2f 毫秒\n", (t_end - t_start) * 1000.0);
    }

    // 6. 汇总所有进程筛选结果
    int* recv_counts = NULL;
    int* displs = NULL;
    size_t* temp_counts = NULL;
    if (rank == 0) {
        recv_counts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        temp_counts = malloc(size * sizeof(size_t));
    }

    // 收集所有进程的计数
    MPI_Gather(&local_count, 1, MPI_UNSIGNED_LONG_LONG, temp_counts, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    // 计算偏移量并转换类型
    if (rank == 0) {
        size_t total_count = 0;
        for (int i = 0; i < size; i++) {
            recv_counts[i] = (int)temp_counts[i];
            displs[i] = (i == 0) ? 0 : displs[i - 1] + recv_counts[i - 1];
            total_count += temp_counts[i];
        }
        printf("最终通过筛选的明密文对数：%zu\n", total_count);
    }

    // 收集所有密文对
    StatePair* all_pairs = NULL;
    if (rank == 0) {
        size_t total_count = 0;
        for (int i = 0; i < size; i++) {
            total_count += temp_counts[i];
        }
        all_pairs = malloc(total_count * sizeof(StatePair));
    }

    MPI_Gatherv(local_pairs, (int)local_count, pair_type, all_pairs, recv_counts, displs, pair_type, 0, MPI_COMM_WORLD);

    // 7. 猜测密钥
    if (rank == 0) {
        size_t total_count = 0;
        for (int i = 0; i < size; i++) {
            total_count += temp_counts[i];
        }

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

        guess_subkeys(matrix_sbox_table, all_pairs, (int)total_count, index_diff_word, &max_val, &max_index_count, max_index);

        // 8. 导出全部候选子密钥
        FILE* fp = fopen("cand_key_0.txt", "w");   // 打开文件，"w" 表示写入模式

        // 将数组元素逐个写入文件
        for (int i = 0; i < max_index_count; i++) {
            fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
        }

        fclose(fp);  // 关闭文件
        printf("数组已成功写入到 cand_key_0.txt\n");
    }

    // 清理内存
    free(local_pairs);
    free(batch);
    if (rank == 0) {
        free(recv_counts);
        free(displs);
        free(temp_counts);
        free(all_pairs);
    }
    MPI_Type_free(&pair_type);
    MPI_Finalize();

    return 0;
}