#define _CRT_SECURE_NO_WARNINGS
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
#define round 100

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
void generate_structure(StatePair* pairs, int num_pairs, const int* diff_word, const int* diff_bits);
void encrypt_and_filter(StatePair* pairs, int num_pairs, const KeySchedule key_schedule, uint8_t*** aes_table, const uint8_t diff_output[16], StatePair* pairs_cand1, StatePair* pairs_cand, int* pairs_cand_count, uint8_t* matrix_sol_table);
void guess_subkeys(const uint8_t matrix_sol_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index);


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

// 生成随机明文结构
void generate_structure(StatePair* pairs, int num_pairs, const int* diff_word,
    const int* diff_bits) {
    for (int i = 0; i < num_pairs; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            uint8_t v = rand() & 0xFF;

            if (diff_word[j] == 1) {
                uint8_t delta = 0;
                for (int bit = 0; bit < 8; bit++) {
                    if (diff_bits[j * 8 + bit]) {
                        delta |= (1 << bit);
                    }
                }
                pairs[i].first[j] = v;
                pairs[i].second[j] = v ^ delta;
            }
            else {
                pairs[i].first[j] = pairs[i].second[j] = v;
            }
        }
    }
}

// 单独将加密过程和筛选密文过程分离出来
void encrypt_and_filter(StatePair* pairs, int num_pairs, const KeySchedule key_schedule, uint8_t*** aes_table, const uint8_t diff_output[16], StatePair* pairs_cand1, StatePair* pairs_cand, int* pairs_cand_count, uint8_t* matrix_sol_table) {
    int max_threads = omp_get_max_threads();
    // 每个线程用私有buffer, 以免冲突
    StatePair** local_cand1_array = malloc(max_threads * sizeof(StatePair*));
    StatePair** local_cand_array = malloc(max_threads * sizeof(StatePair*));
    int* local_count = calloc(max_threads, sizeof(int));

    for (int t = 0; t < max_threads; ++t) {
        local_cand1_array[t] = malloc(num_pairs * sizeof(StatePair));
        local_cand_array[t] = malloc(num_pairs * sizeof(StatePair));
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int local_idx = 0;
        StatePair* local_cand1 = local_cand1_array[tid];
        StatePair* local_cand = local_cand_array[tid];

        #pragma omp for schedule(static)
        for (int i = 0; i < num_pairs; ++i) {
            State c1, c2;
            memcpy(c1, pairs[i].first, STATE_SIZE);
            memcpy(c2, pairs[i].second, STATE_SIZE);

            // XOR with initial key
            for (int j = 0; j < 16; ++j) {
                c1[j] ^= key_schedule[0][j];
                c2[j] ^= key_schedule[0][j];
            }

            // Perform AES rounds
            for (int r = 0; r < round; ++r) {
                State c1_temp, c2_temp;
                memcpy(c1_temp, c1, STATE_SIZE);
                memcpy(c2_temp, c2, STATE_SIZE);

                for (int j = 0; j < 4; ++j) {
                    c1[j] = aes_table[0][j][c1_temp[0]] ^ aes_table[1][j][c1_temp[5]] ^
                        aes_table[2][j][c1_temp[10]] ^ aes_table[3][j][c1_temp[15]];
                    c1[j + 4] = aes_table[0][j][c1_temp[4]] ^ aes_table[1][j][c1_temp[9]] ^
                        aes_table[2][j][c1_temp[14]] ^ aes_table[3][j][c1_temp[3]];
                    c1[j + 8] = aes_table[0][j][c1_temp[8]] ^ aes_table[1][j][c1_temp[13]] ^
                        aes_table[2][j][c1_temp[2]] ^ aes_table[3][j][c1_temp[7]];
                    c1[j + 12] = aes_table[0][j][c1_temp[12]] ^ aes_table[1][j][c1_temp[1]] ^
                        aes_table[2][j][c1_temp[6]] ^ aes_table[3][j][c1_temp[11]];

                    c2[j] = aes_table[0][j][c2_temp[0]] ^ aes_table[1][j][c2_temp[5]] ^
                        aes_table[2][j][c2_temp[10]] ^ aes_table[3][j][c2_temp[15]];
                    c2[j + 4] = aes_table[0][j][c2_temp[4]] ^ aes_table[1][j][c2_temp[9]] ^
                        aes_table[2][j][c2_temp[14]] ^ aes_table[3][j][c2_temp[3]];
                    c2[j + 8] = aes_table[0][j][c2_temp[8]] ^ aes_table[1][j][c2_temp[13]] ^
                        aes_table[2][j][c2_temp[2]] ^ aes_table[3][j][c2_temp[7]];
                    c2[j + 12] = aes_table[0][j][c2_temp[12]] ^ aes_table[1][j][c2_temp[1]] ^
                        aes_table[2][j][c2_temp[6]] ^ aes_table[3][j][c2_temp[11]];
                }

                for (int j = 0; j < 16; ++j) {
                    c1[j] ^= key_schedule[r + 1][j];
                    c2[j] ^= key_schedule[r + 1][j];
                }
            }

            // Final round
            State c1_temp, c2_temp;
            memcpy(c1_temp, c1, STATE_SIZE);
            memcpy(c2_temp, c2, STATE_SIZE);
            sub_bytes(c1, matrix_sol_table);
            shift_rows(c1);
            for (int j = 0; j < 16; ++j) {
                c1[j] ^= key_schedule[round + 1][j];
            }

            sub_bytes(c2, matrix_sol_table);
            shift_rows(c2);
            for (int j = 0; j < 16; ++j) {
                c2[j] ^= key_schedule[round + 1][j];
            }

            // Compare output difference
            bool flag = true;
            for (int j = 0; j < 16; ++j) {
                if ((c1[j] ^ c2[j]) != diff_output[j]) {
                    flag = false;
                    break;
                }
            }

            if (flag) {
                memcpy(local_cand1[local_idx].first, c1_temp, STATE_SIZE);
                memcpy(local_cand1[local_idx].second, c2_temp, STATE_SIZE);
                memcpy(local_cand[local_idx].first, c1, STATE_SIZE);
                memcpy(local_cand[local_idx].second, c2, STATE_SIZE);
                local_idx++;
            }
        }
        local_count[tid] = local_idx;
    }

    // 合并所有线程buffer到全局数组
    int total = 0;
    for (int t = 0; t < max_threads; ++t) {
        if (local_count[t] > 0) {
            memcpy(&pairs_cand1[total], local_cand1_array[t], local_count[t] * sizeof(StatePair));
            memcpy(&pairs_cand[total], local_cand_array[t], local_count[t] * sizeof(StatePair));
            total += local_count[t];
        }
        free(local_cand1_array[t]);
        free(local_cand_array[t]);
    }
    free(local_cand1_array);
    free(local_cand_array);
    free(local_count);

    *pairs_cand_count = total;
}

// 单独将猜测子密钥过程分离出来
void guess_subkeys(const uint8_t matrix_sol_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index) {
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    uint8_t diff_det1 = diff_output[index_inv[index_diff_word[0]]] ^ 1;
    uint8_t diff_det2 = diff_output[index_inv[index_diff_word[1]]] ^ 1;

    int* cand_key_match = (int*)calloc(1 << 16, sizeof(int));

    for (int i = 0; i < (1 << 16); ++i) {
        uint8_t key1 = i & 0xFF;
        uint8_t key2 = (i >> 8) & 0xFF;
        int num = 0;

        for (int j = 0; j < pairs_cand_count; ++j) {
            uint8_t ct1, ct2, ct3, ct4;

            ct1 = pairs_cand[j].first[index_inv[index_diff_word[0]]] ^ key1;
            ct2 = pairs_cand[j].first[index_inv[index_diff_word[1]]] ^ key2;
            ct3 = pairs_cand[j].second[index_inv[index_diff_word[0]]] ^ key1;
            ct4 = pairs_cand[j].second[index_inv[index_diff_word[1]]] ^ key2;

            if (ct1 == 0) ct1 = 1;
            else if (ct1 == 1) ct1 = 0;
            ct1 = matrix_sol_table[ct1];

            if (ct2 == 0) ct2 = 1;
            else if (ct2 == 1) ct2 = 0;
            ct2 = matrix_sol_table[ct2];

            if (ct3 == 0) ct3 = 1;
            else if (ct3 == 1) ct3 = 0;
            ct3 = matrix_sol_table[ct3];

            if (ct4 == 0) ct4 = 1;
            else if (ct4 == 1) ct4 = 0;
            ct4 = matrix_sol_table[ct4];

            bool match = true;
            if ((ct1 ^ ct3) != diff_det1) match = false;
            if ((ct2 ^ ct4) != diff_det2) match = false;

            if (match) num++;
        }
        cand_key_match[i] = num;
    }

    for (int i = 1; i < (1 << 16); ++i) {
        if (cand_key_match[i] > *max_val) {
            *max_val = cand_key_match[i];
        }
    }
    for (int i = 0; i < (1 << 16); ++i) {
        if (cand_key_match[i] == *max_val) {
            max_index[*max_index_count] = i;
            (*max_index_count)++;
        }
    }

    printf("Maximum element: %d\n", *max_val);
    printf("Maximum index count: %d\n", *max_index_count);
    for (int i = 0; i < *max_index_count; ++i) {
        printf("Index of maximum element: %d\n", max_index[i]);

        uint8_t cand_key1 = max_index[i] & 0xFF;
        uint8_t cand_key2 = (max_index[i] >> 8) & 0xFF;
        printf("Candidate key1: 0x%02X, Candidate key2: 0x%02X\n", cand_key1, cand_key2);
    }
}

int main() {
    // 1. 读取必要的数据
    int rows1, cols1, rows2, cols2, rows3, cols3;
    int** diff_sol = read_data_int("diff_sol_0.txt", 10, &rows1, &cols1);
    int** matrix_sol = read_data_int("matrix_sol.txt", 10, &rows2, &cols2);
    int** diff_word_temp = read_data_int("diff_word_0.txt", 10, &rows3, &cols3);
    printf("不等式读取完毕\n");

    // 2. 初始化变量和数据结构
    srand((unsigned int)time(NULL));
    int* diff_word = diff_word_temp[0];
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };

    // 3. 计算 diff_output
    State diff_output;
    compute_diff_output(diff_word, diff_sol, matrix_sol, index, diff_output);
    print_state(diff_output);

    // 4. 构建查找表
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    build_lookup_tables(matrix_sol, matrix_sol_table, xtime_table);

    // 5. 构建AES表
    uint8_t*** aes_table = NULL;
    build_table_aes(matrix_sol_table, xtime_table, &aes_table);

    // 6. 生成主密钥以及密钥拓展
    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    KeySchedule key_schedule;
    key_expansion(key, matrix_sol_table, key_schedule);

    // 7. 生成明文对结构
    StatePair* pairs = (StatePair*)check_malloc(MAX_PAIRS * sizeof(StatePair));
    generate_structure(pairs, MAX_PAIRS, diff_word_temp[0], diff_sol[0]);
    printf("构造明文对数：%d\n", MAX_PAIRS);

    // 8. 加密过程 & 筛选 pairs
    StatePair* pairs_cand = (StatePair*)check_malloc(MAX_PAIRS * sizeof(StatePair));
    StatePair* pairs_cand1 = (StatePair*)check_malloc(MAX_PAIRS * sizeof(StatePair));
    int pairs_cand_count = 0;

    clock_t start = clock();
	encrypt_and_filter(pairs, MAX_PAIRS, key_schedule, aes_table, diff_output, pairs_cand1, pairs_cand, &pairs_cand_count, matrix_sol_table);
    clock_t end = clock();
    double duration = ((double)(end - start)) / CLOCKS_PER_SEC * 1000.0;
    printf("循环耗时: %.2f 毫秒\n", duration);
    printf("满足输出差分的密文对数: %d\n", pairs_cand_count);

    // 9. 猜测密钥
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

	guess_subkeys(matrix_sol_table, pairs_cand, pairs_cand_count, diff_output, index_diff_word, &max_val, &max_index_count, max_index);

    //uint8_t right_key1 = key_schedule[round + 1][index_inv[index_diff_word[0]]];
    //uint8_t right_key2 = key_schedule[round + 1][index_inv[index_diff_word[1]]];
    //printf("Right key1: 0x%02X, Right key2: 0x%02X\n", right_key1, right_key2); // 大写十六进制输出：0xAB
    for (int i = 0; i < round + 2; ++i) {
        for (int j = 0; j < 16; ++j) {
            printf("0x%02X, ", key_schedule[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    //// 10. 导出全部候选子密钥
    //FILE* fp = fopen("cand_key_0.txt", "w");   // 打开文件，"w" 表示写入模式

    //// 将数组元素逐个写入文件
    //for (int i = 0; i < max_index_count; i++) {
    //    fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
    //}

    //fclose(fp);  // 关闭文件
    //printf("数组已成功写入到 cand_key_0.txt\n");

    // 11. 释放内存
    for (int i = 0; i < rows1; i++) free(diff_sol[i]);
    for (int i = 0; i < rows2; i++) free(matrix_sol[i]);
    for (int i = 0; i < rows3; i++) free(diff_word_temp[i]);
    free(diff_sol);
    free(matrix_sol);
    free(diff_word_temp);
    free(pairs);
    free(pairs_cand);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            free(aes_table[i][j]);
        }
        free(aes_table[i]);
    }
    free(aes_table);

    return 0;
}