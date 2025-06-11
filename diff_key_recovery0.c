#define _CRT_SECURE_NO_WARNINGS
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
#define MAX_LINE_LENGTH 1024
#define MAX_PAIRS (1ULL<<32)
#define MAX_FILTERED_PAIRS 10000  // 最多在内存中保留多少对通过筛选的明密文对
#define round 100

// Type definitions
typedef uint8_t State[STATE_SIZE];
typedef uint8_t KeySchedule[round + 2][16];

// Structure definitions
typedef struct {
    State first;
    State second;
} StatePair;

typedef struct {
    uint32_t state;
    int rank;
} RandomState;

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
int** read_data_int(const char* filename, int base, int* rows, int* cols);
void init_random_state(RandomState* rs, int rank);
uint32_t next_random(RandomState* rs);
void compute_diff_output(int* diff_word, int** diff_sol, int** matrix_sol, const int index[], uint8_t diff_output[16]);
void print_state(State diff);
void build_lookup_tables(int** matrix_sol, uint8_t matrix_sol_table[256], uint8_t xtime_table[256]);
void sub_bytes(State s, uint8_t* table);
void shift_rows(State state);
void build_table_aes(uint8_t* sbox_table, uint8_t* time_table, uint8_t**** aes_table);
void generate_rcon(uint8_t* rcon);
void key_expansion(const uint8_t* key, uint8_t* table, KeySchedule round_keys);
void generate_pair(StatePair* pair, const int* diff_word, const int* diff_bits, RandomState* rs);
void encrypt(State plaintext, State ciphertext, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table);
bool pass_filter(State c1, State c2, const uint8_t diff_output[16]);
void generate_encrypt_and_filter_stream(size_t total_pairs, const int* diff_word, const int* diff_bits, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table, const uint8_t diff_output[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, RandomState* rs);
void guess_subkeys(const uint8_t matrix_sol_table[256], const StatePair* pairs_cand, int pairs_cand_count, const uint8_t diff_output[16], const int index_diff_word[2], int* max_val, int* max_index_count, int* max_index);
void broadcast_matrix(int** matrix, int rows, int cols, int root, MPI_Comm comm);


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

void init_random_state(RandomState* rs, int rank) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    // 使用纳秒级时间戳提高初始随机性
    rs->state = (uint32_t)(ts.tv_nsec + rank * 1299709);
    rs->rank = rank;
}

uint32_t next_random(RandomState* rs) {
    // 使用线性同余法生成随机数
    rs->state = rs->state * 1103515245 + 12345;
    return rs->state;
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
void generate_pair(StatePair* pair, const int* diff_word, const int* diff_bits, RandomState* rs) {
    for (int j = 0; j < STATE_SIZE; j++) {
        // 使用next_random生成随机数，确保每次调用都能得到不同的随机数
        uint8_t v = next_random(rs) & 0xFF;  // 使用新的随机数生成器
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
void generate_encrypt_and_filter_stream(size_t total_pairs, const int* diff_word, const int* diff_bits, const KeySchedule key_schedule, uint8_t*** aes_table, uint8_t* matrix_sol_table, const uint8_t diff_output[16], StatePair* filtered_pairs, size_t* filtered_count, size_t max_filtered_count, RandomState* rs)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 计算每个进程的工作量
    size_t pairs_per_proc = total_pairs / size;
    size_t start_pair = rank * pairs_per_proc;
    size_t end_pair = (rank == size - 1) ? total_pairs : (rank + 1) * pairs_per_proc;

    // 本地缓冲区 - 使用check_malloc替代malloc
    const size_t BATCH_SIZE = 1024;
    StatePair* local_buf = check_malloc(max_filtered_count * sizeof(StatePair));
    size_t local_count = 0;

    //// 进程特定的随机数种子
    //uint32_t seed = (uint32_t)(time(NULL) ^ (rank * 7919));

    // 临时数组
    StatePair batch[BATCH_SIZE];
    State c1[BATCH_SIZE], c2[BATCH_SIZE];

    // 主处理循环
    for (size_t i = start_pair; i < end_pair; i += BATCH_SIZE) {
        size_t batch_end = i + BATCH_SIZE;
        if (batch_end > end_pair) batch_end = end_pair;
        size_t batch_size = batch_end - i;

        // 批量处理
        for (size_t j = 0; j < batch_size; ++j) {
            //generate_pair(&batch[j], diff_word, diff_bits, &seed);
            // 使用新的随机数生成器
            generate_pair(&batch[j], diff_word, diff_bits, rs);  // 修改这行

            encrypt(batch[j].first, c1[j], key_schedule, aes_table, matrix_sol_table);
            encrypt(batch[j].second, c2[j], key_schedule, aes_table, matrix_sol_table);

            if (pass_filter(c1[j], c2[j], diff_output)) {
                if (local_count < max_filtered_count) {
                    memcpy(local_buf[local_count].first, c1[j], STATE_SIZE);
                    memcpy(local_buf[local_count].second, c2[j], STATE_SIZE);
                    local_count++;
                }
            }
        }
    }

    // 创建MPI数据类型用于StatePair
    MPI_Datatype mpi_statepair;
    MPI_Type_contiguous(sizeof(StatePair), MPI_BYTE, &mpi_statepair);
    MPI_Type_commit(&mpi_statepair);

    // 收集所有进程的结果数量
    size_t* all_counts = NULL;
    size_t* displs = NULL;
    if (rank == 0) {
        all_counts = check_malloc(size * sizeof(size_t));  // 使用check_malloc
        displs = check_malloc(size * sizeof(size_t));      // 使用check_malloc
    }

    // 收集每个进程的局部结果数量
    MPI_Gather(&local_count, 1, MPI_UNSIGNED_LONG_LONG, all_counts, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    // 计算位移和总数
    if (rank == 0) {
        displs[0] = 0;
        size_t total = all_counts[0];
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + all_counts[i - 1];
            total += all_counts[i];
        }
        *filtered_count = (total > max_filtered_count) ? max_filtered_count : total;
    }

    // 收集所有结果
    MPI_Gatherv(local_buf, local_count, mpi_statepair, filtered_pairs, all_counts, displs, mpi_statepair, 0, MPI_COMM_WORLD);

    // 清理
    MPI_Type_free(&mpi_statepair);
    free(local_buf);
    if (rank == 0) {
        free(all_counts);
        free(displs);
    }
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

int main(int argc, char** argv) {
    int rank, size;
    double start_time, end_time;

    // 1. 初始化MPI
    MPI_CHECK(MPI_Init(&argc, &argv));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));

    // 记录开始时间
    start_time = MPI_Wtime();

    // 初始化随机数生成器
    RandomState rs;
    init_random_state(&rs, rank);

    if (rank == 0) {
        printf("Running with %d processes\n", size);
        printf("Program initialization completed, starting main computation...\n\n");
    }

    // 这些变量所有进程都需要，所以在外面声明
    // 2. 声明所需变量
    int rows1, cols1, rows2, cols2, rows3, cols3;
    int** diff_sol = NULL, ** matrix_sol = NULL, ** diff_word_temp = NULL;
    int* diff_word = NULL;
    State diff_output;
    uint8_t matrix_sol_table[256] = { 0 };
    uint8_t xtime_table[256] = { 0 };
    uint8_t*** aes_table = NULL;
    KeySchedule key_schedule;
    uint8_t key[] = { 0x3F, 0xA9, 0x72, 0x5C, 0x01, 0xBD, 0xE2, 0x9F, 0x56, 0x11, 0x3A, 0xC4, 0xD8, 0x77, 0x99, 0xAB };
    int index[] = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
    int index_inv[] = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };

    // 3. 主进程读取数据并进行初始计算
    if (rank == 0) {
        // 读取数据
        diff_sol = read_data_int("diff_sol_0.txt", 10, &rows1, &cols1);
        matrix_sol = read_data_int("matrix_sol.txt", 10, &rows2, &cols2);
        diff_word_temp = read_data_int("diff_word_0.txt", 10, &rows3, &cols3);
        diff_word = diff_word_temp[0];

        printf("Data reading completed.\n");

        // 初始计算
        compute_diff_output(diff_word, diff_sol, matrix_sol, index, diff_output);
        printf("Differential output computed:\n");
        print_state(diff_output);

        // 构建查找表
        build_lookup_tables(matrix_sol, matrix_sol_table, xtime_table);
        printf("Lookup tables built.\n");

        // 构建AES表
        build_table_aes(matrix_sol_table, xtime_table, &aes_table);
        printf("AES tables built.\n");

        // 密钥扩展
        key_expansion(key, matrix_sol_table, key_schedule);
        printf("Key schedule completed.\n");

        // 测试随机数生成器
        printf("Testing random number generator:\n");
        for (int i = 0; i < 5; i++) {
            printf("Random %d: %u\n", i, next_random(&rs));
        }
    }

    // 4. 广播必要数据给所有进程
    // 广播维度信息
    MPI_CHECK(MPI_Bcast(&rows1, 1, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(&cols1, 1, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(&rows2, 1, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(&cols2, 1, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(&rows3, 1, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(&cols3, 1, MPI_INT, 0, MPI_COMM_WORLD));

    // 广播矩阵数据
    broadcast_matrix(diff_sol, rows1, cols1, 0, MPI_COMM_WORLD);
    broadcast_matrix(matrix_sol, rows2, cols2, 0, MPI_COMM_WORLD);
    broadcast_matrix(diff_word_temp, rows3, cols3, 0, MPI_COMM_WORLD);

    // 广播计算结果
    MPI_CHECK(MPI_Bcast(diff_output, STATE_SIZE, MPI_UINT8_T, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(matrix_sol_table, 256, MPI_UINT8_T, 0, MPI_COMM_WORLD));
    MPI_CHECK(MPI_Bcast(xtime_table, 256, MPI_UINT8_T, 0, MPI_COMM_WORLD));

    // 广播AES表（需要特殊处理三维数组）
    if (rank != 0) {
        aes_table = (uint8_t***)check_malloc(4 * sizeof(uint8_t**));
        for (int i = 0; i < 4; i++) {
            aes_table[i] = (uint8_t**)check_malloc(4 * sizeof(uint8_t*));
            for (int j = 0; j < 4; j++) {
                aes_table[i][j] = (uint8_t*)check_malloc(256 * sizeof(uint8_t));
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            MPI_CHECK(MPI_Bcast(aes_table[i][j], 256, MPI_UINT8_T, 0, MPI_COMM_WORLD));
        }
    }

    // 广播密钥扩展结果
    MPI_CHECK(MPI_Bcast(key_schedule, sizeof(KeySchedule), MPI_BYTE, 0, MPI_COMM_WORLD));

    if (rank == 0) {
        printf("All data broadcast completed. Proceeding to main computation...\n");
    }

    // 5. 同步所有进程
    MPI_Barrier(MPI_COMM_WORLD);

    // 6. 生成明文 & 加密过程 & 筛选 pairs
    StatePair* filtered_pairs = NULL;
    size_t filtered_count = 0;
    if (rank == 0) {
        filtered_pairs = malloc(MAX_FILTERED_PAIRS * sizeof(StatePair));
    }

    // 在generate_encrypt_and_filter_stream函数中使用随机数生成器
    // 需要修改函数签名以接收RandomState参数, 添加随机数生成器参数&rs
    generate_encrypt_and_filter_stream(MAX_PAIRS, diff_word, diff_sol[0], key_schedule, aes_table, matrix_sol_table, diff_output, filtered_pairs, &filtered_count, MAX_FILTERED_PAIRS, &rs);

    // 7. 只在主进程处理结果和写文件
    if (rank == 0) {
        printf("最终通过筛选的明密文对数：%zu\n", filtered_count);
        /*clock_t end = clock();
        double duration = ((double)(end - start)) / CLOCKS_PER_SEC * 1000.0;
        printf("循环耗时: %.2f 毫秒\n", duration);*/
        // ... 结果处理和文件写入 ...

        // 猜测密钥
        int index_diff_word[2];
        int index_diff_word_count = 0;
        for (int i = 0; i < 16; ++i) {
            if (diff_word[i]) {
                index_diff_word[index_diff_word_count++] = i;
            }
        }
        int max_val = 0;
        int* max_index = check_malloc((1 << 16) * sizeof(int));  // 使用check_malloc
        int max_index_count = 0;

        guess_subkeys(matrix_sol_table, filtered_pairs, filtered_count, diff_output, index_diff_word, &max_val, &max_index_count, max_index);

        // 导出全部候选子密钥
        FILE* fp = fopen("cand_key_0.txt", "w");   // 打开文件，"w" 表示写入模式

        // 将数组元素逐个写入文件
        for (int i = 0; i < max_index_count; i++) {
            fprintf(fp, "%d\n", max_index[i]);  // 每个元素占一行
        }

        fclose(fp);  // 关闭文件
        printf("数组已成功写入到 cand_key_0.txt\n");

        free(max_index);  // 确实需要这行代码来释放内存

        end_time = MPI_Wtime();
        printf("总运行时间: %.2f 秒\n", end_time - start_time);
    }
    
    // 8. 清理内存
    if (rank == 0) {
        // 释放内存
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
    }
    else {
        // 非主进程也需要释放aes_table
        if (aes_table != NULL) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    free(aes_table[i][j]);
                }
                free(aes_table[i]);
            }
            free(aes_table);
        }
    }

    MPI_CHECK(MPI_Finalize());
    return 0;
}