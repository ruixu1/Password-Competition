#ifndef XOSHIRO256PLUSPLUS_H
#define XOSHIRO256PLUSPLUS_H

#include <stdint.h>

typedef struct {
    uint64_t s[4];
} xoshiro256plusplus_state;

// 初始化，推荐每个线程用不同seed和线程号
void xoshiro256plusplus_init(xoshiro256plusplus_state* state, uint64_t seed, uint64_t thread_id);

// 生成一个64位随机数
uint64_t xoshiro256plusplus_next(xoshiro256plusplus_state* state);

#endif