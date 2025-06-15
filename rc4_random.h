#ifndef RC4_RANDOM_H
#define RC4_RANDOM_H

#include <stdint.h>

typedef struct {
    uint8_t S[256];
    uint8_t a, b;
} RC4State;

// 初始化RC4状态，seed建议用足够随机的32位数
void rc4_init(RC4State* state, uint32_t seed);

// 输出下一个32位随机数
uint32_t rc4_next(RC4State* state);

#endif