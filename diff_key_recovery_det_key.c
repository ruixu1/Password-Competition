#define _CRT_SECURE_NO_WARNINGS
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
#define round 100
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
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
int** read_data_int(const char* filename, int base, int* rows, int* cols);
void compute_diff_output(int* diff_word, int** diff_sol, int** matrix_sol, const int index[], uint8_t diff_output[16]);
void print_state(State diff);
void build_lookup_tables(int** matrix_sol, uint8_t matrix_sol_table[256], uint8_t xtime_table[256]);
void sub_bytes(State s, uint8_t* table);
void shift_rows(State state);
void build_table_aes(uint8_t* sbox_table, uint8_t* time_table, uint8_t**** aes_table);
void generate_rcon(uint8_t* rcon);
void key_expansion(const uint8_t* key, uint8_t* table, KeySchedule round_keys);
void extract_index_diff_word(int* diff_word, int index_diff_word[2], int* index_diff_word_count);
void aes_key_schedule_invert(const State last_round_key, uint8_t* table, KeySchedule round_keys);
void encrypt(State plaintext, State ciphertext, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table);



// Function declarations
void* check_malloc(size_t size) {
    void* ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    return ptr;
}

// File reading functions
int** read_data_int(const char* filename, int base, int* rows, int* cols) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "No such file: %s\n", filename);
        exit(1);
    }

    char line[MAX_LINE_LENGTH];
    int max_cols = 0;
    *rows = 0;

    // Count rows and maximum columns
    while (fgets(line, sizeof(line), file)) {
        int col_count = 0;
        char* token = strtok(line, " \t\n");
        while (token) {
            col_count++;
            token = strtok(NULL, " \t\n");
        }
        if (col_count > max_cols) max_cols = col_count;
        (*rows)++;
    }
    *cols = max_cols;

    // Allocate memory for the matrix
    int** matrix = (int**)check_malloc(*rows * sizeof(int*));
    for (int i = 0; i < *rows; i++) {
        matrix[i] = (int*)check_malloc(max_cols * sizeof(int));
    }

    // Reset file pointer and read data
    rewind(file);
    int row = 0;
    while (fgets(line, sizeof(line), file)) {
        int col = 0;
        char* token = strtok(line, " \t\n");
        while (token) {
            matrix[row][col] = (int)strtol(token, NULL, base);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }

    fclose(file);
    return matrix;
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

// 给定diff_word_0，生成对应索引集合diff_word_temp_0
void extract_index_diff_word(const State diff_word, int index_diff_word[2], int* index_diff_word_count) {
    *index_diff_word_count = 0;
    for (int i = 0; i < 16; ++i) {
        if (diff_word[i]) {
            index_diff_word[*index_diff_word_count] = i;
            (*index_diff_word_count)++;
        }
    }
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

// 逆密钥编排，需要S盒查找表和Rcon表
void aes_key_schedule_invert(const State last_round_key, const uint8_t matrix_sbox_table[256], KeySchedule round_keys) {
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

    // 初始化最后一轮密钥
    memcpy(round_keys[round + 1], last_round_key, STATE_SIZE);

    for (int r = round + 1; r > 0; --r) {
        // 逆推每一轮
        uint8_t* current = round_keys[r];
        uint8_t* prev = round_keys[r - 1];

        // 逆推最后一列（第13-15字节）
        for (int j = 15; j >= 4; --j) {
            prev[j] = current[j] ^ current[j - 4];
        }

        // 逆推前4字节
        uint8_t temp[4];
        temp[0] = prev[12];
        temp[1] = prev[13];
        temp[2] = prev[14];
        temp[3] = prev[15];

        // RotWord
        uint8_t t = temp[0];
        temp[0] = temp[1]; temp[1] = temp[2]; temp[2] = temp[3]; temp[3] = t;

        // SubWord
        temp[0] = matrix_sbox_table[temp[0]];
        temp[1] = matrix_sbox_table[temp[1]];
        temp[2] = matrix_sbox_table[temp[2]];
        temp[3] = matrix_sbox_table[temp[3]];

        // Rcon
        temp[0] ^= rcon[r];

        prev[0] = current[0] ^ temp[0];
        prev[1] = current[1] ^ temp[1];
        prev[2] = current[2] ^ temp[2];
        prev[3] = current[3] ^ temp[3];
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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint64_t total = ((uint64_t)1 << 32);
    uint64_t chunk = (total + size - 1) / size;
    uint64_t my_start = rank * chunk;
    uint64_t my_end = my_start + chunk;
    if (my_end > total) my_end = total;

    // 1. 读取必要的数据
    int rows[8], cols[8];
    int** cand_key_0 = read_data_int("cand_key_0.txt", 10, &rows[0], &cols[0]);
    int** cand_key_1 = read_data_int("cand_key_1.txt", 10, &rows[1], &cols[1]);
    int** cand_key_2 = read_data_int("cand_key_2.txt", 10, &rows[2], &cols[2]);
    int** cand_key_3 = read_data_int("cand_key_3.txt", 10, &rows[3], &cols[3]);
    int** cand_key_4 = read_data_int("cand_key_4.txt", 10, &rows[4], &cols[4]);
    int** cand_key_5 = read_data_int("cand_key_5.txt", 10, &rows[5], &cols[5]);
    int** cand_key_6 = read_data_int("cand_key_6.txt", 10, &rows[6], &cols[6]);
    int** cand_key_7 = read_data_int("cand_key_7.txt", 10, &rows[7], &cols[7]);
    printf("不等式读取完毕\n");

    // 2. 初始化变量和数据结构
    srand((unsigned int)time(NULL));
    int i, j;
    State diff_word_0 = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    State diff_word_1 = { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    State diff_word_2 = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
    State diff_word_3 = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1 };
    State diff_word_4 = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
    State diff_word_5 = { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 };
    State diff_word_6 = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 };
    State diff_word_7 = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 };
    uint8_t matrix_sol[] = { 0x97, 0x92, 0x94, 0x98, 0x10, 0x20, 0x40, 0x80 };
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    int index_diff_word[8][2];
    int index_diff_word_count = 0;
    uint64_t global_seed = 0xC7E6A3B2F19D4E58ULL; // 全局种子，用于xoshiro256plusplus初始化

    extract_index_diff_word(diff_word_0, index_diff_word[0], &index_diff_word_count);
    extract_index_diff_word(diff_word_1, index_diff_word[1], &index_diff_word_count);
    extract_index_diff_word(diff_word_2, index_diff_word[2], &index_diff_word_count);
    extract_index_diff_word(diff_word_3, index_diff_word[3], &index_diff_word_count);
    extract_index_diff_word(diff_word_4, index_diff_word[4], &index_diff_word_count);
    extract_index_diff_word(diff_word_5, index_diff_word[5], &index_diff_word_count);
    extract_index_diff_word(diff_word_6, index_diff_word[6], &index_diff_word_count);
    extract_index_diff_word(diff_word_7, index_diff_word[7], &index_diff_word_count);

    // 3. 构建查找表
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t matrix_sbox_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, matrix_sbox_table, xtime_table);

    uint8_t aes_table[16][256] = { { 0x00 } };
    build_table_aes(matrix_sbox_table, xtime_table, aes_table);

    // 5. 生成主密钥以及密钥拓展
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sbox_table, rcon, key_schedule);

    // 6. 随机生成10组明文并且加密得到对应密文
	uint8_t plain[10][STATE_SIZE];
    uint8_t cipher[10][STATE_SIZE];
    for (i = 0; i < 10; ++i) {
        xoshiro256plusplus_state prng;
        xoshiro256plusplus_init(&prng, global_seed, i + 1);
        uint64_t r1 = xoshiro256plusplus_next(&prng);
        uint64_t r2 = xoshiro256plusplus_next(&prng);

        for (j = 0; j < 8; ++j) {
            plain[i][j] = (r1 >> (8 * j)) & 0xFF;
            plain[i][j + 8] = (r2 >> (8 * j)) & 0xFF;
        }

		encrypt(plain[i], cipher[i], key_schedule, aes_table, matrix_sbox_table, index);
    }

    // 7. 得到2^{32}种全部可能的候选子密钥，并且和生成的明文加密验证得到正确子密钥
    double t1 = MPI_Wtime();
    int found = 0;
    uint64_t result_index = 0;
    uint8_t result_key[16] = { 0 };
    State key_cand_right = { 0xE3, 0x6C, 0xDB, 0xBF, 0x7F, 0x05, 0x91, 0x05, 0xCF, 0xCC, 0x94, 0x27, 0x02, 0x3D, 0xB9, 0xC2 };

    for (uint64_t i = my_start; i < my_end; ++i) {
        // ... 组装key_cand_local，逆推密钥，验证10组明密文 ...
        State key_cand_local;  // 每个线程自己的key_cand
        KeySchedule key_schedule_inv;
        uint8_t cipher_cand[10][STATE_SIZE];

        int index[8];
        for (int j = 0; j < 8; ++j) {
            index[j] = (i >> (4 * j)) & 0xF; // 注意加括号
        }
        uint8_t temp[16];
        temp[0] = cand_key_0[index[0]][0] & 0xFF;
        temp[1] = (cand_key_0[index[0]][0] >> 8) & 0xFF;
        temp[2] = cand_key_1[index[1]][0] & 0xFF;
        temp[3] = (cand_key_1[index[1]][0] >> 8) & 0xFF;
        temp[4] = cand_key_2[index[2]][0] & 0xFF;
        temp[5] = (cand_key_2[index[2]][0] >> 8) & 0xFF;
        temp[6] = cand_key_3[index[3]][0] & 0xFF;
        temp[7] = (cand_key_3[index[3]][0] >> 8) & 0xFF;
        temp[8] = cand_key_4[index[4]][0] & 0xFF;
        temp[9] = (cand_key_4[index[4]][0] >> 8) & 0xFF;
        temp[10] = cand_key_5[index[5]][0] & 0xFF;
        temp[11] = (cand_key_5[index[5]][0] >> 8) & 0xFF;
        temp[12] = cand_key_6[index[6]][0] & 0xFF;
        temp[13] = (cand_key_6[index[6]][0] >> 8) & 0xFF;
        temp[14] = cand_key_7[index[7]][0] & 0xFF;
        temp[15] = (cand_key_7[index[7]][0] >> 8) & 0xFF;
        for (int j = 0; j < 8; ++j) {
            for (int k = 0; k < 2; ++k) {
                key_cand_local[index_inv[index_diff_word[j][k]]] = temp[2 * j + k];
            }
        }

        aes_key_schedule_invert(key_cand_local, matrix_sbox_table, key_schedule_inv);

        // 使用key_schedule_inv验证明密文对
        bool match = true;
        for (int j = 0; j < 10; ++j) {
            encrypt(plain[j], cipher_cand[j], key_schedule_inv, aes_table, matrix_sbox_table, index);
            match = true;
            for (int k = 0; k < 16; ++k) {
                if (cipher_cand[j][k] != cipher[j][k]) {
                    match = false;
                    break;
                }
            }
            if (!match) break;
        }

        if (match) {
            // 抢答写全局结果
            {
                if (!found) {
                    found = 1;
                    result_index = i;
                    memcpy(result_key, key_cand_local, 16);
                }
            }
        }
        if (match) {
            found = 1;
            result_index = i;
            memcpy(result_key, key_cand_local, 16);
            break; // 本进程提前退出
        }
    }

    double t2 = MPI_Wtime();

    if (found) {
        printf("Rank %d found key! index=%" PRIu64 "\n", rank, result_index);
        printf("Key: ");
        for (j = 0; j < 16; ++j) printf("0x%02X, ", result_key[j]);
        printf("\n");
        printf("耗时: %.2f ms\n", (t2 - t1) * 1000.0);
    }
    else {
        printf("No matching key found.\n");
    }

    // 8. 释放内存
    for (int i = 0; i < rows[0]; i++) free(cand_key_0[i]);
    for (int i = 0; i < rows[1]; i++) free(cand_key_1[i]);
    for (int i = 0; i < rows[2]; i++) free(cand_key_2[i]);
    for (int i = 0; i < rows[3]; i++) free(cand_key_3[i]);
    for (int i = 0; i < rows[4]; i++) free(cand_key_4[i]);
    for (int i = 0; i < rows[5]; i++) free(cand_key_5[i]);
    for (int i = 0; i < rows[6]; i++) free(cand_key_6[i]);
    for (int i = 0; i < rows[7]; i++) free(cand_key_7[i]);
    free(cand_key_0);
    free(cand_key_1);
    free(cand_key_2);
    free(cand_key_3);
    free(cand_key_4);
    free(cand_key_5);
    free(cand_key_6);
    free(cand_key_7);

    MPI_Finalize();
    return 0;

    return 0;
}