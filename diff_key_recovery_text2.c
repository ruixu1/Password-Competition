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
//#define MAX_PAIRS (1ULL<<32)
#define MAX_PAIRS 200000
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
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
void generate_pair(StatePair* pair, const int* diff_word, const int* diff_bits, uint32_t* seed);
void encrypt(State plaintext, State ciphertext, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table);
bool pass_filter(State c1, State c2, const uint8_t diff_output[16]);
void generate_encrypt_and_filter_stream(size_t num_pairs, const int* diff_word, const int* diff_bits, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table, const uint8_t diff_output[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count);
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

// 添加在函数声明之后，main函数之前
static inline uint32_t fast_rand(uint32_t* state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

// 生成随机明文结构
void generate_pair(StatePair* pair, const int* diff_word, const int* diff_bits, uint32_t* seed) {  // 添加seed参数
    for (int j = 0; j < STATE_SIZE; j++) {
        uint8_t v = fast_rand(seed) & 0xFF;  // 使用fast_rand替换rand
		// 使用fast_rand生成随机数，确保每次调用都能得到不同的随机数
        if (diff_word[j] == 1) {
            uint8_t delta = 0;
            for (int bit = 0; bit < 8; bit++) {
                if (diff_bits[j * 8 + bit]) {
                    delta |= (1 << bit);
                }
            }
            pair->first[j] = v;
            pair->second[j] = v ^ delta;
        }
        else {
            pair->first[j] = pair->second[j] = v;
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

// 判断一对加密后是否满足差分（你原本的逻辑）
bool pass_filter(State c1, State c2, const uint8_t diff_output[16]) {
    for (int i = 0; i < STATE_SIZE; ++i) {
        if ((c1[i] ^ c2[i]) != diff_output[i]) return false;
    }
    return true;
}

// 单独将加密过程和筛选密文过程分离出来
void generate_encrypt_and_filter_stream(size_t total_pairs, const int* diff_word, const int* diff_bits, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table, const uint8_t diff_output[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count) {
    // 限制线程数为物理核心数
    int max_threads = omp_get_max_threads();
    if (max_threads > 32) max_threads = 32;  // 调整为你的服务器物理核心数

    // 使用更大的线程本地buffer
    size_t thread_bufsize = (max_filtered_count / max_threads) + 1000;
    StatePair** thread_buffers = malloc(max_threads * sizeof(StatePair*));
    size_t* thread_counts = calloc(max_threads, sizeof(size_t));

    for (int t = 0; t < max_threads; ++t) {
        thread_buffers[t] = malloc(thread_bufsize * sizeof(StatePair));
    }

    const size_t BATCH_SIZE = 1024;  // 批处理大小

    #pragma omp parallel num_threads(max_threads)
    {
        int tid = omp_get_thread_num();
        StatePair* local_buf = thread_buffers[tid];
        size_t local_count = 0;

        // 线程独立的随机数种子
        uint32_t seed = (uint32_t)(time(NULL) ^ (tid * 7919));

        // 批处理临时数组
        StatePair batch[BATCH_SIZE];
        State c1[BATCH_SIZE], c2[BATCH_SIZE];

        #pragma omp for schedule(dynamic, 1024)
        for (size_t i = 0; i < total_pairs; i += BATCH_SIZE) {
            size_t batch_end = i + BATCH_SIZE;
            if (batch_end > total_pairs) batch_end = total_pairs;

            // 1. 批量生成明文对
            for (size_t j = 0; j < batch_end - i; ++j) {
                generate_pair(&batch[j], diff_word, diff_bits, &seed);
            }

            // 2. 批量加密
            for (size_t j = 0; j < batch_end - i; ++j) {
                encrypt(batch[j].first, c1[j], key_schedule, aes_table, matrix_sol_table);
                encrypt(batch[j].second, c2[j], key_schedule, aes_table, matrix_sol_table);
            }

            // 3. 批量过滤
            for (size_t j = 0; j < batch_end - i; ++j) {
                if (pass_filter(c1[j], c2[j], diff_output)) {
                    if (local_count < thread_bufsize) {
                        local_buf[local_count++] = batch[j];
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
void guess_subkeys(const uint8_t matrix_sol_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index) {
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };
    uint8_t diff_det1 = diff_output[index_inv[index_diff_word[0]]] ^ 1;
    uint8_t diff_det2 = diff_output[index_inv[index_diff_word[1]]] ^ 1;

    printf("match diff1: 0x%02X, match diff2: 0x%02X\n", diff_det1, diff_det2);
    printf("\n");

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
    /*for (int i = 0; i < *max_index_count; ++i) {
        printf("Index of maximum element: %d\n", max_index[i]);

        uint8_t cand_key1 = max_index[i] & 0xFF;
        uint8_t cand_key2 = (max_index[i] >> 8) & 0xFF;
        printf("Candidate key1: 0x%02X, Candidate key2: 0x%02X\n", cand_key1, cand_key2);
    }*/
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
    printf("diff_output : ");
    print_state(diff_output);
    printf("\n");

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

    // 7. 生成明文 & 加密过程 & 筛选 pairs
    StatePair* filtered_pairs = malloc(MAX_FILTERED_PAIRS * sizeof(StatePair));
    size_t filtered_count = 0;

    clock_t start = clock();
    generate_encrypt_and_filter_stream(
        MAX_PAIRS, diff_word, diff_sol[0], key_schedule,
        aes_table, matrix_sol_table, diff_output,
        filtered_pairs, &filtered_count, MAX_FILTERED_PAIRS
    );
    printf("最终通过筛选的明密文对数：%zu\n", filtered_count);
    clock_t end = clock();
    double duration = ((double)(end - start)) / CLOCKS_PER_SEC * 1000.0;
    printf("循环耗时: %.2f 毫秒\n", duration);
    for (int i = 0; i < filtered_count; ++i) {
        for (int j = 0; j < 16; ++j) {
            uint8_t temp = filtered_pairs[i].first[j] ^ key_schedule[round + 1][j];
            printf("%02X ", temp);
        }
        printf("\n");
        for (int j = 0; j < 16; ++j) {
            uint8_t temp = filtered_pairs[i].second[j] ^ key_schedule[round + 1][j];
            printf("%02X ", temp);
        }
        printf("\n");
    }

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

	guess_subkeys(matrix_sol_table, filtered_pairs, filtered_count, diff_output, index_diff_word, &max_val, &max_index_count, max_index);

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
    for (int i = 0; i < max_index_count; i++) {
        fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
    }

    fclose(fp);  // 关闭文件
    printf("数组已成功写入到 cand_key_0.txt\n");

    // 11. 释放内存
    for (int i = 0; i < rows1; i++) free(diff_sol[i]);
    for (int i = 0; i < rows2; i++) free(matrix_sol[i]);
    for (int i = 0; i < rows3; i++) free(diff_word_temp[i]);
    free(diff_sol);
    free(matrix_sol);
    free(diff_word_temp);
    free(filtered_pairs);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            free(aes_table[i][j]);
        }
        free(aes_table[i]);
    }
    free(aes_table);

    return 0;
}