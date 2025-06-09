#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <omp.h>

// Constants
#define STATE_SIZE 16   // 128位
#define BYTE_SIZE 8
#define MAX_LINE_LENGTH 1024
#define MAX_PAIRS 200000
#define round 20

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

// 拆出 diff_output 计算函数
void compute_diff_output(int* diff_word, int** diff_sol, int** matrix_sol, const int index[], uint8_t diff_output[16]) {
    uint8_t diff_temp1[128] = { 0 };
    uint8_t diff_temp2[128] = { 0 };

    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 8; ++j) {
            diff_temp1[i * 8 + j] = diff_sol[0][i * 8 + j];
        }
        if (diff_word[i]) {
            diff_temp1[i * 8] ^= 1;
        }
    }

    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 8; ++j) {
            for (int k = 0; k < 8; ++k) {
                diff_temp2[i * 8 + j] ^= matrix_sol[j][k] * diff_temp1[i * 8 + k];
            }
        }
    }

    for (int i = 0; i < 16; ++i) {
        uint8_t byte = 0;
        for (int j = 0; j < 8; j++) {
            byte |= (diff_temp2[index[i] * 8 + j] & 0x1) << j;
        }
        diff_output[i] = byte;
    }
}

void print_state(State diff) {
    for (int i = 0; i < 16; ++i) {
        printf("%02X ", diff[i]);
    }
    printf("\n");
}

// 拆出 lookup tables 构建
void build_lookup_tables(int** matrix_sol, uint8_t matrix_sol_table[256], uint8_t xtime_table[256]) {
    for (int i = 0; i < 256; i++) {
        int temp[8] = { 0 }, temp1[8] = { 0 };
        for (int j = 0; j < 8; j++) {
            temp[j] = (i >> j) & 1;
        }
        for (int j = 0; j < 8; j++) {
            temp1[j] = 0;
            for (int k = 0; k < 8; k++) {
                temp1[j] ^= matrix_sol[j][k] * temp[k];
            }
        }
        uint8_t temp_value = 0;
        for (int j = 0; j < 8; j++) {
            temp_value |= (temp1[j] << j);
        }
        matrix_sol_table[i] = temp_value;
    }

    for (int i = 0; i < 256; i++) {
        xtime_table[i] = (i << 1) ^ ((i & 0x80) ? 0x1B : 0x00);
    }
}

// AES操作相关函数
void sub_bytes(State s, uint8_t* table) {
    for (int i = 0; i < STATE_SIZE; i++) {
        if (s[i] == 0) s[i] = 1;
        else if (s[i] == 1) s[i] = 0;
        s[i] = table[s[i]];
    }
}

void shift_rows(State state) {
    State temp;
    memcpy(temp, state, STATE_SIZE);

    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            int src_idx = r + 4 * ((c + r) % 4);
            state[r + 4 * c] = temp[src_idx];
        }
    }
}

// 表格构建函数
void build_table_aes(uint8_t* sbox_table, uint8_t* time_table, uint8_t**** aes_table) {
    // 分配内存
    *aes_table = (uint8_t***)check_malloc(4 * sizeof(uint8_t**));
    for (int i = 0; i < 4; i++) {
        (*aes_table)[i] = (uint8_t**)check_malloc(4 * sizeof(uint8_t*));
        for (int j = 0; j < 4; j++) {
            (*aes_table)[i][j] = (uint8_t*)check_malloc(256 * sizeof(uint8_t));
        }
    }

    for (int i = 0; i < 256; i++) {
        // S盒变换
        uint8_t x = i;
        if (x == 0) x = 1;
        else if (x == 1) x = 0;
        x = sbox_table[x];

        // 列混淆和行移位
        (*aes_table)[0][0][i] = time_table[x];
        (*aes_table)[0][1][i] = x;
        (*aes_table)[0][2][i] = x;
        (*aes_table)[0][3][i] = time_table[x] ^ x;

        (*aes_table)[1][0][i] = time_table[x] ^ x;
        (*aes_table)[1][1][i] = time_table[x];
        (*aes_table)[1][2][i] = x;
        (*aes_table)[1][3][i] = x;

        (*aes_table)[2][0][i] = x;
        (*aes_table)[2][1][i] = time_table[x] ^ x;
        (*aes_table)[2][2][i] = time_table[x];
        (*aes_table)[2][3][i] = x;

        (*aes_table)[3][0][i] = x;
        (*aes_table)[3][1][i] = x;
        (*aes_table)[3][2][i] = time_table[x] ^ x;
        (*aes_table)[3][3][i] = time_table[x];
    }
}

// Rcon生成函数
void generate_rcon(uint8_t* rcon) {
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
void key_expansion(const uint8_t* key, uint8_t* table, KeySchedule round_keys) {
    uint8_t rcon[round + 2];
    generate_rcon(rcon);

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
        for (int j = 0; j < 4; ++j) {
            if (temp[j] == 0) temp[j] = 1;
            else if (temp[j] == 1) temp[j] = 0;
            temp[j] = table[temp[j]];
        }

        // Add Rcon
        temp[0] ^= rcon[r];

        // Generate first 4 bytes
        for (int j = 0; j < 4; ++j) {
            round_keys[r][j] = round_keys[r - 1][j] ^ temp[j];
        }

        // Generate remaining 12 bytes
        for (int j = 4; j < 16; ++j) {
            round_keys[r][j] = round_keys[r - 1][j] ^ round_keys[r][j - 4];
        }
    }
}

// 给定diff_word_0，生成对应索引集合diff_word_temp_0
void extract_index_diff_word(int* diff_word, int index_diff_word[2], int* index_diff_word_count) {
    *index_diff_word_count = 0;
    for (int i = 0; i < 16; ++i) {
        if (diff_word[i]) {
            index_diff_word[*index_diff_word_count] = i;
            (*index_diff_word_count)++;
        }
    }
}

// 需要S盒查找表和Rcon表
void aes_key_schedule_invert(const State last_round_key, uint8_t* table, KeySchedule round_keys) {
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
        for (int j = 0; j < 4; ++j) temp[j] = prev[j + 12];

        // RotWord
        uint8_t t = temp[0];
        temp[0] = temp[1]; temp[1] = temp[2]; temp[2] = temp[3]; temp[3] = t;

        // SubWord
        for (int j = 0; j < 4; ++j) {
            if (temp[j] == 0) temp[j] = 1;
            else if (temp[j] == 1) temp[j] = 0;
            temp[j] = table[temp[j]];
        }

        // Rcon
        temp[0] ^= rcon[r];

        for (int j = 0; j < 4; ++j) {
            prev[j] = current[j] ^ temp[j];
        }
    }
}

// round轮aes加密
void encrypt(State plaintext, State ciphertext, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table) {
    memcpy(ciphertext, plaintext, STATE_SIZE);
    // XOR with initial key
    for (int j = 0; j < 16; ++j) {
        ciphertext[j] ^= key_schedule[0][j];
    }

    // Perform AES rounds
    for (int r = 0; r < round; ++r) {
        State c_temp;
        memcpy(c_temp, ciphertext, STATE_SIZE);

        for (int j = 0; j < 4; ++j) {
            ciphertext[j] = aes_table[0][j][c_temp[0]] ^ aes_table[1][j][c_temp[5]] ^
                aes_table[2][j][c_temp[10]] ^ aes_table[3][j][c_temp[15]];
            ciphertext[j + 4] = aes_table[0][j][c_temp[4]] ^ aes_table[1][j][c_temp[9]] ^
                aes_table[2][j][c_temp[14]] ^ aes_table[3][j][c_temp[3]];
            ciphertext[j + 8] = aes_table[0][j][c_temp[8]] ^ aes_table[1][j][c_temp[13]] ^
                aes_table[2][j][c_temp[2]] ^ aes_table[3][j][c_temp[7]];
            ciphertext[j + 12] = aes_table[0][j][c_temp[12]] ^ aes_table[1][j][c_temp[1]] ^
                aes_table[2][j][c_temp[6]] ^ aes_table[3][j][c_temp[11]];
        }

        for (int j = 0; j < 16; ++j) {
            ciphertext[j] ^= key_schedule[r + 1][j];
        }
    }

    // Final round
    sub_bytes(ciphertext, matrix_sol_table);
    shift_rows(ciphertext);
    for (int j = 0; j < 16; ++j) {
        ciphertext[j] ^= key_schedule[round + 1][j];
    }
}


int main() {
    // 1. 读取必要的数据
    int rows[17], cols[17];
    int** matrix_sol = read_data_int("matrix_sol.txt", 10, &rows[0], &cols[0]);
    int** diff_word_temp_0 = read_data_int("diff_word_0.txt", 10, &rows[1], &cols[1]);
    int** diff_word_temp_1 = read_data_int("diff_word_1.txt", 10, &rows[2], &cols[2]);
    int** diff_word_temp_2 = read_data_int("diff_word_2.txt", 10, &rows[3], &cols[3]);
    int** diff_word_temp_3 = read_data_int("diff_word_3.txt", 10, &rows[4], &cols[4]);
    int** diff_word_temp_4 = read_data_int("diff_word_4.txt", 10, &rows[5], &cols[5]);
    int** diff_word_temp_5 = read_data_int("diff_word_5.txt", 10, &rows[6], &cols[6]);
    int** diff_word_temp_6 = read_data_int("diff_word_6.txt", 10, &rows[7], &cols[7]);
    int** diff_word_temp_7 = read_data_int("diff_word_7.txt", 10, &rows[8], &cols[8]);
    int** cand_key_0 = read_data_int("cand_key_0.txt", 10, &rows[9], &cols[9]);
    int** cand_key_1 = read_data_int("cand_key_1.txt", 10, &rows[10], &cols[10]);
    int** cand_key_2 = read_data_int("cand_key_2.txt", 10, &rows[11], &cols[11]);
    int** cand_key_3 = read_data_int("cand_key_3.txt", 10, &rows[12], &cols[12]);
    int** cand_key_4 = read_data_int("cand_key_4.txt", 10, &rows[13], &cols[13]);
    int** cand_key_5 = read_data_int("cand_key_5.txt", 10, &rows[14], &cols[14]);
    int** cand_key_6 = read_data_int("cand_key_6.txt", 10, &rows[15], &cols[15]);
    int** cand_key_7 = read_data_int("cand_key_7.txt", 10, &rows[16], &cols[16]);
    printf("不等式读取完毕\n");

    // 2. 初始化变量和数据结构
    srand((unsigned int)time(NULL));
    int* diff_word_0 = diff_word_temp_0[0];
    int* diff_word_1 = diff_word_temp_1[0];
    int* diff_word_2 = diff_word_temp_2[0];
    int* diff_word_3 = diff_word_temp_3[0];
    int* diff_word_4 = diff_word_temp_4[0];
    int* diff_word_5 = diff_word_temp_5[0];
    int* diff_word_6 = diff_word_temp_6[0];
    int* diff_word_7 = diff_word_temp_7[0];
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    int index_diff_word[8][2];
    int index_diff_word_count = 0;

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
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, xtime_table);

    // 4. 构建AES表
    uint8_t*** aes_table = NULL;
    build_table_aes(matrix_sol_table, xtime_table, &aes_table);

    // 5. 生成主密钥以及密钥拓展
    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sol_table, key_schedule);

    // 6. 随机生成10组明文并且加密得到对应密文
	uint8_t plain[10][STATE_SIZE];
    uint8_t cipher[10][STATE_SIZE];
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 16; j++) {
            plain[i][j] = rand() & 0xFF;
        }
		encrypt(plain[i], cipher[i], key_schedule, aes_table, matrix_sol_table);
    }

    // 7. 得到2^{32}种全部可能的候选子密钥， 并且和生成的明文加密验证得到正确子密钥
     // 7. 得到2^{32}种全部可能的候选子密钥， 并且和生成的明文加密验证得到正确子密钥
    clock_t start = clock();
    State key_cand_right = { 0xE3, 0x6C, 0xDB, 0xBF, 0x7F, 0x05, 0x91, 0x05, 0xCF, 0xCC, 0x94, 0x27, 0x02, 0x3D, 0xB9, 0xC2 };
    State key_cand;

    for (uint64_t i = 0; i < ((uint64_t)1 << 32); ++i) {
        int index[8];
        for (int j = 0; j < 8; ++j) {
            index[j] = (i >> 4 * j) & 0xF; // 取低4位作为索引
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
                key_cand[index_inv[index_diff_word[j][k]]] = temp[2 * j + k];
            }
        }

        KeySchedule key_schedule_inv;
        aes_key_schedule_invert(key_cand, matrix_sol_table, key_schedule_inv);

        // 使用key_schedule_inv验证明密文对
        uint8_t cipher_cand[10][STATE_SIZE];
        bool match = true;
        for (int j = 0; j < 10; ++j) {
            encrypt(plain[j], cipher_cand[j], key_schedule_inv, aes_table, matrix_sol_table);
            match = true;
            for (int k = 0; k < 16; ++k) {
                if (cipher_cand[j][k] != cipher[j][k]) {
                    match = false;
                    break;
                }
            }
            if (!match) {
                break;
            }
        }

        if (match) {
            printf("The correct index is : %" PRIu64 "\n", i);
            printf("The correct sub key is : ");
            for (int j = 0; j < 16; ++j) {
                printf("0x%02X, ", key_cand[j]);
            }
            break;
        }
    }

    clock_t end = clock();
    double duration = ((double)(end - start)) / CLOCKS_PER_SEC * 1000.0;
    printf("循环耗时: %.2f 毫秒\n", duration);

    // 8. 释放内存
    for (int i = 0; i < rows[0]; i++) free(matrix_sol[i]);
    for (int i = 0; i < rows[1]; i++) free(diff_word_temp_0[i]);
    for (int i = 0; i < rows[2]; i++) free(diff_word_temp_1[i]);
    for (int i = 0; i < rows[3]; i++) free(diff_word_temp_2[i]);
    for (int i = 0; i < rows[4]; i++) free(diff_word_temp_3[i]);
    for (int i = 0; i < rows[5]; i++) free(diff_word_temp_4[i]);
    for (int i = 0; i < rows[6]; i++) free(diff_word_temp_5[i]);
    for (int i = 0; i < rows[7]; i++) free(diff_word_temp_6[i]);
    for (int i = 0; i < rows[8]; i++) free(diff_word_temp_7[i]);
    for (int i = 0; i < rows[9]; i++) free(cand_key_0[i]);
    for (int i = 0; i < rows[10]; i++) free(cand_key_1[i]);
    for (int i = 0; i < rows[11]; i++) free(cand_key_2[i]);
    for (int i = 0; i < rows[12]; i++) free(cand_key_3[i]);
    for (int i = 0; i < rows[13]; i++) free(cand_key_4[i]);
    for (int i = 0; i < rows[14]; i++) free(cand_key_5[i]);
    for (int i = 0; i < rows[15]; i++) free(cand_key_6[i]);
    for (int i = 0; i < rows[16]; i++) free(cand_key_7[i]);
    free(matrix_sol);
    free(diff_word_temp_0);
    free(diff_word_temp_1);
    free(diff_word_temp_2);
    free(diff_word_temp_3);
    free(diff_word_temp_4);
    free(diff_word_temp_5);
    free(diff_word_temp_6);
    free(diff_word_temp_7);
    free(cand_key_0);
    free(cand_key_1);
    free(cand_key_2);
    free(cand_key_3);
    free(cand_key_4);
    free(cand_key_5);
    free(cand_key_6);
    free(cand_key_7);
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            free(aes_table[i][j]);
        }
        free(aes_table[i]);
    }
    free(aes_table);

    return 0;
}