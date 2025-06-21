#define _CRT_SECURE_NO_WARNINGS
#include "xoshiro256plusplus.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
//#include <omp.h>
#include <mpi.h>

// Constants
#define STATE_SIZE 16   // 128位
#define BYTE_SIZE 8
#define BATCH_SIZE 4096  // 将批处理大小定义为常量
#define MAX_LINE_LENGTH 1024
#define MAX_PAIRS (1ULL<<32)
//#define MAX_PAIRS 10000000
#define round 101
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
#define AES_TABLE_ROW(OUT_IDX, T0, T1, T2, T3, IN_IDX_0, IN_IDX_1, IN_IDX_2, IN_IDX_3) \
    ciphertext[OUT_IDX] = aes_table[T0][c_temp[IN_IDX_0]] ^ \
                          aes_table[T1][c_temp[IN_IDX_1]] ^ \
                          aes_table[T2][c_temp[IN_IDX_2]] ^ \
                          aes_table[T3][c_temp[IN_IDX_3]]

// Type definitions
typedef uint8_t State[STATE_SIZE];
typedef uint8_t KeySchedule[round + 3][16];

// Structure definitions
typedef struct {
    State first;
    State second;
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
void* check_malloc(size_t size);
void print_state(State diff);
void broadcast_matrix(int** matrix, int rows, int cols, int root, MPI_Comm comm);
void build_lookup_tables(const uint8_t matrix_sol[8], uint8_t matrix_sol_table[256], uint8_t matrix_sbox_table[256], uint8_t xtime_table[256]);
void compute_diff_output(const uint8_t diff_word[16], const uint8_t diff_sol[16], const uint8_t matrix_sol_table[256], const int index[16], State diff_output);
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
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint8_t diff_bites[16], const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const uint8_t diff_output[16], const int index[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed, int mpi_rank, int mpi_size);
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

// 广播矩阵数据到所有进程
void broadcast_matrix(int** matrix, int rows, int cols, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    // 广播矩阵维度
    int dims[2] = { rows, cols };
    MPI_CHECK(MPI_Bcast(dims, 2, MPI_INT, root, comm));

    // 非root进程分配内存
    if (rank != root) {
        matrix = (int**)check_malloc(dims[0] * sizeof(int*));
        for (int i = 0; i < dims[0]; i++) {
            matrix[i] = (int*)check_malloc(dims[1] * sizeof(int));
        }
    }

    // 广播数据
    for (int i = 0; i < dims[0]; i++) {
        MPI_CHECK(MPI_Bcast(matrix[i], dims[1], MPI_INT, root, comm));
    }
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
void compute_diff_output(const uint8_t diff_word[16], const uint8_t diff_sol[16], const uint8_t matrix_sol_table[256], const int index[16], State diff_output) {
    for (int i = 0; i < 16; ++i) {
        uint8_t diff_temp = diff_sol[index[i]] ^ diff_word[index[i]];
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
void key_expansion(const uint8_t key[16], const uint8_t matrix_sbox_table[256], const uint8_t rcon[round + 3], KeySchedule round_keys) {
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
        ciphertext[j] ^= key_schedule[round + 2][j];
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
void generate_encrypt_and_filter_stream(size_t total_pairs, const uint8_t diff_bites[16], const KeySchedule key_schedule, const uint8_t aes_table[16][256], const uint8_t matrix_sbox_table[256], const uint8_t diff_output[16], const int index[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, uint64_t global_seed, int mpi_rank, int mpi_size) {
    size_t i, j;
    // 1. 计算本进程需要处理的pair区间
    size_t pairs_per_rank = (total_pairs + mpi_size - 1) / mpi_size;
    size_t my_start = mpi_rank * pairs_per_rank;
    size_t my_end = my_start + pairs_per_rank;
    if (my_end > total_pairs) my_end = total_pairs;

    // 2. 分配对齐的本地缓冲区
    StatePair* batch = (StatePair*)aligned_alloc(64, BATCH_SIZE * sizeof(StatePair));
    State* c1 = (State*)aligned_alloc(64, BATCH_SIZE * sizeof(State));
    State* c2 = (State*)aligned_alloc(64, BATCH_SIZE * sizeof(State));

    // 用于记录本地筛选结果的计数器
    size_t local_count = 0;

    // 3. 初始化PRNG
    xoshiro256plusplus_state prng;
    xoshiro256plusplus_init(&prng, global_seed, mpi_rank);

    // 4. 主循环：批量生成、加密、筛选
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
                if (local_count < max_filtered_count) {  // 使用传入的max_filtered_count
                    memcpy(filtered_pairs[local_count].first, c1[j], STATE_SIZE);
                    memcpy(filtered_pairs[local_count].second, c2[j], STATE_SIZE);
                    local_count++;
                }
            }
        }
    }

    // 5. 返回筛选结果数量
    *filtered_count = local_count;

    free(batch);
    free(c1);
    free(c2);
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
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    State diff_sol = { 0x00, 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4f, 0x00, 0x00, 0x00, 0x00 };
    State diff_plain = { 0x01, 0x00, 0x9f, 0x00, 0x00, 0x26, 0x00, 0xf7, 0x01, 0x00, 0x9f, 0x00, 0x00, 0x26, 0x00, 0xf7 };
    State diff_word = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    uint8_t matrix_sol[] = { 0x97, 0x92, 0x94, 0x98, 0x10, 0x20, 0x40, 0x80 };
    uint64_t global_seed = 0xC7E6A3B2F19D4E58ULL; // 全局种子，用于xoshiro256plusplus初始化

    // 2. 构建查找表
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t matrix_sbox_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, matrix_sbox_table, xtime_table);

    uint8_t aes_table[16][256] = { { 0x00 } };
    build_table_aes(matrix_sbox_table, xtime_table, aes_table);

    // 3. 计算 diff_output
    State diff_output;
    compute_diff_output(diff_word, diff_sol, matrix_sol_table, index, diff_output);
    print_state(diff_output);

    // 4. 生成主密钥以及密钥拓展
    // 生成轮常数rcon
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sbox_table, rcon, key_schedule);

    // 5. 各进程独立生成、加密、筛选
    // 创建MPI派生类型
    MPI_Datatype pair_type;
    MPI_Type_contiguous(32, MPI_BYTE, &pair_type);
    MPI_Type_commit(&pair_type);

    MPI_Barrier(MPI_COMM_WORLD); // 所有进程同步，开始计时
    double t_start = MPI_Wtime();

    StatePair* local_pairs = malloc(MAX_FILTERED_PAIRS * sizeof(StatePair));
    size_t local_count = 0;

    generate_encrypt_and_filter_stream(MAX_PAIRS, diff_plain, key_schedule, aes_table, matrix_sbox_table, diff_output, index, local_pairs, &local_count, MAX_FILTERED_PAIRS, global_seed, rank, size);

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

        guess_subkeys(matrix_sbox_table, all_pairs, (int)total_count, diff_output, index_diff_word, &max_val, &max_index_count, max_index);

        // 8. 导出全部候选子密钥
        FILE* fp = fopen("cand_key_1.txt", "w");   // 打开文件，"w" 表示写入模式

        // 将数组元素逐个写入文件
        for (int i = 0; i < max_index_count; i++) {
            fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
        }

        fclose(fp);  // 关闭文件
        printf("数组已成功写入到 cand_key_1.txt\n");
    }

    // 清理内存
    free(local_pairs);
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