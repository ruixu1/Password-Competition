#define _CRT_SECURE_NO_WARNINGS
#include <bits/stdc++.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include<stdio.h>
#include <iomanip>
#include <sstream>
#include <set>
#include <array>
#include <bitset>
#include <random>
#include <ctime>
#include <chrono>
#include <iterator>


using namespace std;

// AES参数
const int STATE_SIZE = 16;   // 128位
const int BYTE_SIZE = 8;

// 类型定义
using State = array<uint8_t, STATE_SIZE>;
using Bit128 = bitset<128>;
using KeySchedule = vector<vector<uint8_t>>;
#define uint32 unsigned int


// 用来读取文件
vector<vector<uint32>> ReadData(string filename, int base) {
	//  Can read integers with sign with various length, base can take 10 or 16
	ifstream in(filename);
	if (!in) cout << "No such file!" << endl;
	string s;
	vector<vector<string>> res_str;
	int i, j;
	while (getline(in, s)) {
		string t;
		//  Filter data to store integers (include hexadecimal code), space and x (for 16)
		for (i = 0; i < s.size(); i++) {
			if ((s[i] >= 48 && s[i] <= 57) || s[i] == 32 || s[i] == 120 || s[i] >= 97 && s[i] <= 102 || s[i] == 45) t += s[i];
		}

		j = 0;
		vector<string> tmp;
		for (i = 0; i < t.size(); i++) {
			if (t[i] == 32) {   //  deal with space
				string tmp2 = t.substr(j, i - j + 1);
				tmp.push_back(tmp2);
				j = i + 1;
			}
		}

		if (t[t.size() - 1] != 32) {
			string tmp2 = t.substr(j);
			tmp.push_back(tmp2);
		}
		res_str.push_back(tmp);
	}

	vector<vector<uint32>> result;
	for (i = 0; i < res_str.size(); i++) {
		vector<uint32> tmp;
		for (j = 0; j < res_str[i].size(); j++) {
			const char* nptr = res_str[i][j].c_str();
			char* endptr = NULL;
			errno = 0;
			int val = strtol(nptr, &endptr, base);
			tmp.push_back(val);
		}
		result.push_back(tmp);

	}

	return result;
}


vector<vector<int>> ReadDataInt(string filename, int base) {
	//  Can read integers with sign with various length, base can take 10 or 16
	ifstream in(filename);
	if (!in) cout << "No such file!" << endl;
	string s;
	vector<vector<string>> res_str;
	int i, j;
	while (getline(in, s)) {
		string t;
		//  Filter data to store integers (include hexadecimal code), space and x (for 16)
		for (i = 0; i < s.size(); i++) {
			if ((s[i] >= 48 && s[i] <= 57) || s[i] == 32 || s[i] == 120 || s[i] >= 97 && s[i] <= 102 || s[i] == 45) t += s[i];
		}

		j = 0;
		vector<string> tmp;
		for (i = 0; i < t.size(); i++) {
			if (t[i] == 32) {   //  deal with space
				string tmp2 = t.substr(j, i - j + 1);
				tmp.push_back(tmp2);
				j = i + 1;
			}
		}

		if (t[t.size() - 1] != 32) {
			string tmp2 = t.substr(j);
			tmp.push_back(tmp2);
		}
		res_str.push_back(tmp);
	}

	vector<vector<int>> result;
	for (i = 0; i < res_str.size(); i++) {
		vector<int> tmp;
		for (j = 0; j < res_str[i].size(); j++) {
			const char* nptr = res_str[i][j].c_str();
			char* endptr = NULL;
			errno = 0;
			int val = strtol(nptr, &endptr, base);
			tmp.push_back(val);
		}
		result.push_back(tmp);

	}

	return result;
}


void SubBytes(State& s, vector<uint8_t> table) {
	for (auto& b : s) {
		if (b == 0) b = 1;
		else if (b == 1) b = 0;

		b = table[b];
	}
}


// 行移位 ShiftRows
void ShiftRows(State& state) {
	State temp = state;
	for (int r = 0; r < 4; ++r) {
		for (int c = 0; c < 4; ++c) {
			state[r + 4 * c] = temp[r + 4 * ((c + r) % 4)];
		}
	}
}


// 建立四个 8 进 32 出的表格，用于快速实现aes的一列变换
vector<vector<vector<uint8_t>>> build_table_aes(vector<uint8_t> sbox_table, vector<uint8_t> time_table) {
	vector<vector<vector<uint8_t>>> aes_table(4, vector<vector<uint8_t>>(4, vector<uint8_t>(256, 0)));
	for (int i = 0; i < 256; ++i) {
		// 首先经过S盒变换
		uint8_t x = i;
		if (x == 0) x = 1;
		if (x == 1) x = 0;
		x = sbox_table[x];
		
		// 然后进行列混淆，行移位自己手动区分
		aes_table[0][0][i] = time_table[x];
		aes_table[0][1][i] = x;
		aes_table[0][2][i] = x;
		aes_table[0][3][i] = time_table[x] ^ x;

		aes_table[1][0][i] = time_table[x] ^ x;
		aes_table[1][1][i] = time_table[x];
		aes_table[1][2][i] = x;
		aes_table[1][3][i] = x;

		aes_table[2][0][i] = x;
		aes_table[2][1][i] = time_table[x] ^ x;
		aes_table[2][2][i] = time_table[x];
		aes_table[2][3][i] = x;

		aes_table[3][0][i] = x;
		aes_table[3][1][i] = x;
		aes_table[3][2][i] = time_table[x] ^ x;
		aes_table[3][3][i] = time_table[x];
	}
	return aes_table;
}


// 生成 Rcon 常量
vector<uint8_t> generate_rcon(int round_count) {
	vector<uint8_t> rcon(round_count + 1);
	rcon[0] = 0x00;  // 第0个是占位
	uint8_t value = 0x01;
	for (int i = 1; i <= round_count; ++i) {
		rcon[i] = value;
		value <<= 1;
		if (value & 0x100) {
			value ^= 0x11B;  // 多项式模 GF(2^8)
		}
	}
	return rcon;
}


// 密钥扩展主函数：输入16字节主密钥，输出round_count轮子密钥
KeySchedule KeyExpansion(const vector<uint8_t>& key, int round_count, vector<uint8_t> table) {
	KeySchedule round_keys(round_count + 1, vector<uint8_t>(16));
	vector<uint8_t> rcon = generate_rcon(round_count);

	// 拷贝初始密钥作为第0轮密钥
	for (int i = 0; i < 16; ++i) {
		round_keys[0][i] = key[i];
	}

	for (int round = 1; round <= round_count; ++round) {
		vector<uint8_t> temp(4);

		// 获取上一轮最后4字节
		for (int j = 0; j < 4; ++j) {
			temp[j] = round_keys[round - 1][12 + j];
		}

		// RotWord: 左循环移位
		rotate(temp.begin(), temp.begin() + 1, temp.end());

		// SubWord: 应用 S-box
		for (int j = 0; j < 4; ++j) {
			if (temp[j] == 0) temp[j] = 1;
			else if (temp[j] == 1) temp[j] = 0;

			temp[j] = table[temp[j]];
		}

		// 添加 Rcon 常量（只作用于第一个字节）
		temp[0] ^= rcon[round];

		// 生成前4字节
		for (int j = 0; j < 4; ++j) {
			round_keys[round][j] = round_keys[round - 1][j] ^ temp[j];
		}

		// 后12字节
		for (int j = 4; j < 16; ++j) {
			round_keys[round][j] = round_keys[round - 1][j] ^ round_keys[round][j - 4];
		}
	}

	return round_keys;
}


// 随机明文结构生成：构造满足 diff_word 中定义的字节差分模式的明文对
vector<pair<State, State>> generate_structure(int num_pairs, const vector<int>& diff_word, const vector<int>& diff_bits) {
	random_device rd;                   // 用于获取硬件随机种子
	mt19937 gen(rd());                  // 构造一个高质量伪随机数生成器
	uniform_int_distribution<> dis(0, 255);  // 均匀分布：生成 0 到 255 之间的随机数

	vector<pair<State, State>> pairs;   // 存放明文对 (P1, P2)

	for (int i = 0; i < num_pairs; ++i) {
		State a{}, b{};                 // 初始化两个 16 字节明文状态，默认为全 0

		// 遍历每个字节位置，根据 diff_word 决定是否加差分
		for (int j = 0; j < STATE_SIZE; ++j) {
			uint8_t v = dis(gen);   // 随机生成一个字节

			if (diff_word[j] == 1) {
				// 如果该位置是活跃的（差分参与传递）
				// 从 bit-level 差分中提取当前字节的差分掩码
				uint8_t delta = 0;
				for (int bit = 0; bit < 8; ++bit) {
					if (diff_bits[j * 8 + bit]) {
						delta |= (1 << bit);
					}
				}
				a[j] = v;   // 明文1的这个字节
				b[j] = v ^ delta;  // 按真实差分掩码设置 b[j]
			}
			else {
				// 非活跃字节保持相同（不参与差分）
				a[j] = b[j] = v;
			}
		}

		// 将构造好的明文对加入结果列表中
		pairs.emplace_back(a, b);
	}

	return pairs;  // 返回构造好的明文对结构
}


// 替代原来的函数：将 AES 状态转为 bitset<128>
bitset<128> state_to_bitset(const State& s) {
	bitset<128> b;
	for (int i = 0; i < STATE_SIZE; ++i) {
		for (int j = 0; j < BYTE_SIZE; ++j) {
			b[i * 8 + j] = (s[i] >> j) & 1;
		}
	}
	return b;
}


int main() {
	// 其中diff_sol是比特差分模式，matrix_sol是得到的仿射变换矩阵, diff_word_temp是字节级别差分模式，matrix_aes是AES的线性层128*128矩阵
	vector<vector<int>> diff_sol = ReadDataInt("diff_sol_0.txt", 10);
	vector<vector<int>> matrix_sol = ReadDataInt("matrix_sol.txt", 10);
	vector<vector<int>> diff_word_temp = ReadDataInt("diff_word_0.txt", 10);
	cout << "不等式读取完毕" << endl;

	vector<int> diff_word = diff_word_temp[0];  // 取第一行作为差分模式
	
	// 代表行移位的索引映射
	vector<int> index = { 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 };
	// 代表逆行移位的索引映射
	vector<int> index_inv = { 0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3 };

	// 计算目标密文差分, 即区分器往后拓展一轮后的输出差分
	vector<int> diff_temp1(128, 0), diff_temp2(128, 0);
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 8; ++j) {
			diff_temp1[i * 8 + j] = diff_sol[0][i * 8 + j];
		}
		if (diff_word[i]) {
			diff_temp1[i * 8] = diff_sol[0][i * 8] ^ 1;
		}
	}
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 8; ++j) {
			for (int k = 0; k < 8; ++k) {
				diff_temp2[i * 8 + j] ^= matrix_sol[j][k] * diff_temp1[i * 8 + k];
			}
		}
	}

	bitset<128> diff_sol_bitset, diff_output;
	for (int i = 0; i < 128; ++i) {
		diff_sol_bitset[i] = diff_sol[0][i];
	}
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 8; ++j) {
			diff_output[i * 8 + j] = diff_temp2[index[i] * 8 + j];
		}
	}

	// 仿射变换矩阵matrix_sol的查找表
	vector<uint8_t> matrix_sol_table;
	vector<int> temp(8), temp1(8);
	for (int i = 0; i < 256; ++i) {
		for (int j = 0; j < 8; ++j) {
			temp[j] = (i >> j) & 1;
		}
		for (int j = 0; j < 8; ++j) {
			temp1[j] = 0;
			for (int k = 0; k < 8; ++k) {
				temp1[j] ^= matrix_sol[j][k] * temp[k];
			}
		}

		// 将 temp1 看作一个二进制数
		uint8_t temp_value = 0;
		for (int j = 0; j < 8; ++j) {
			temp_value |= (temp1[j] << j);
		}
		matrix_sol_table.push_back(temp_value);
	}

	// 构建2*x的查找表
	vector<uint8_t> xtime_table(256);
	for (int i = 0; i < 256; ++i) {
		xtime_table[i] = (i << 1) ^ ((i & 0x80) ? 0x1B : 0x00);
	}
	vector<vector<vector<uint8_t>>> aes_table = build_table_aes(matrix_sol_table, xtime_table);

	// 随机主密钥（此版本加密函数未使用，可扩展）
	vector<uint8_t> key(16);
	random_device rd;
	for (auto& b : key) b = rd() & 0xFF;

	int round = 4 * 5; // AES 轮数（此处为简化版本，实际 AES 为 10 轮）
	
	// 扩展密钥（此版本加密函数未使用，可扩展）
	KeySchedule key_schedule = KeyExpansion(key, round + 1, matrix_sol_table);

	// 生成结构
	auto pairs = generate_structure(200000, diff_word, diff_sol[0]);  // 可调数量
	cout << "构造明文对数：" << pairs.size() << endl;

	int match_count = 0, t = 0;

	// 开始计时
	auto start = chrono::high_resolution_clock::now();

	vector<pair<State, State>> pairs_cand;

	for (auto& [p1, p2] : pairs) {
		// 模拟4轮AES（只做SubBytes）
		State c1 = p1, c2 = p2;
		// 首先异或主密钥
		for (int i = 0; i < 16; ++i) {
			c1[i] ^= key_schedule[0][i];
			c2[i] ^= key_schedule[0][i];
		}

		// 进行 round 轮的 AES 加密
		for (int r = 0; r < round; ++r) {
			State c1_temp = c1, c2_temp = c2;
			for (int i = 0; i < 4; ++i) {
				c1[i] = aes_table[0][i][c1_temp[0]] ^ aes_table[1][i][c1_temp[5]] ^ aes_table[2][i][c1_temp[10]] ^ aes_table[3][i][c1_temp[15]];
				c1[i + 4] = aes_table[0][i][c1_temp[4]] ^ aes_table[1][i][c1_temp[9]] ^ aes_table[2][i][c1_temp[14]] ^ aes_table[3][i][c1_temp[3]];
				c1[i + 8] = aes_table[0][i][c1_temp[8]] ^ aes_table[1][i][c1_temp[13]] ^ aes_table[2][i][c1_temp[2]] ^ aes_table[3][i][c1_temp[7]];
				c1[i + 12] = aes_table[0][i][c1_temp[12]] ^ aes_table[1][i][c1_temp[1]] ^ aes_table[2][i][c1_temp[6]] ^ aes_table[3][i][c1_temp[11]];
				c2[i] = aes_table[0][i][c2_temp[0]] ^ aes_table[1][i][c2_temp[5]] ^ aes_table[2][i][c2_temp[10]] ^ aes_table[3][i][c2_temp[15]];
				c2[i + 4] = aes_table[0][i][c2_temp[4]] ^ aes_table[1][i][c2_temp[9]] ^ aes_table[2][i][c2_temp[14]] ^ aes_table[3][i][c2_temp[3]];
				c2[i + 8] = aes_table[0][i][c2_temp[8]] ^ aes_table[1][i][c2_temp[13]] ^ aes_table[2][i][c2_temp[2]] ^ aes_table[3][i][c2_temp[7]];
				c2[i + 12] = aes_table[0][i][c2_temp[12]] ^ aes_table[1][i][c2_temp[1]] ^ aes_table[2][i][c2_temp[6]] ^ aes_table[3][i][c2_temp[11]];
			}

			for (int i = 0; i < 16; ++i) {
				c1[i] ^= key_schedule[r + 1][i];
			}

			for (int i = 0; i < 16; ++i) {
				c2[i] ^= key_schedule[r + 1][i];
			}
		}

		
		// 进行最后一轮（无列混淆操作）， SubBytes 和 ShiftRows
		SubBytes(c1, matrix_sol_table);
		ShiftRows(c1);
		for (int i = 0; i < 16; ++i) {
			c1[i] ^= key_schedule[round + 1][i];
		}

		SubBytes(c2, matrix_sol_table);
		ShiftRows(c2);
		for (int i = 0; i < 16; ++i) {
			c2[i] ^= key_schedule[round + 1][i];
		}

		// 转换为比特，比较是否匹配 diff_output
		// 转换为比特，比较是否匹配 diff_sol[0]
		bitset<128> b1 = state_to_bitset(c1);
		bitset<128> b2 = state_to_bitset(c2);
		bitset<128> d = b1 ^ b2;
		bool flag = (d == diff_output);

		if (flag) {
			pairs_cand.emplace_back(c1, c2);
			match_count++;
		}
	}

	// 结束计时
	auto end = chrono::high_resolution_clock::now();

	// 计算耗时（以毫秒为单位）
	chrono::duration<double, milli> duration = end - start;
	cout << "循环耗时: " << duration.count() << " 毫秒" << endl;

	cout << "满足输出差分的密文对数: " << match_count << endl;

	vector<int> index_diff_word;
	for (int i = 0; i < 16; ++i) {
		if (diff_word[i]) {
			index_diff_word.push_back(i);
		}
	}

	// 遍历两个子密钥所有可能的组合，分别统计每种子密钥的匹配数量
	vector<int> cand_key_match((1 << 16), 0);
	for (int i = 0; i < (1 << 16); ++i) {
		uint8_t key1, key2;
		key1 = i & 0xFF;  // 取低8位
		key2 = (i >> 8) & 0xFF;  // 取高8位
		int num = 0;
		for (auto& [c1, c2] : pairs_cand) {
			uint8_t ct1, ct2, ct3, ct4;
			//异或密钥
			ct1 = c1[index_diff_word[0]] ^ key1;
			ct2 = c1[index_diff_word[1]] ^ key2;
			ct3 = c2[index_diff_word[0]] ^ key1;
			ct4 = c2[index_diff_word[1]] ^ key2;

			//逆S盒
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

			// 检查是否满足差分
			bool match = true;
			uint8_t c_temp;
			c_temp = ct1 ^ ct3;
			if (c_temp != 0x4F) match = false;
			c_temp = ct2 ^ ct4;
			if (c_temp != 0x4F) match = false;

			if (match) num++;
		}
		cand_key_match[i] = num;
	}

	// 1. 使用 max_element 找到最大元素的迭代器
	auto max_it = max_element(cand_key_match.begin(), cand_key_match.end());

	// 2. 使用 distance 计算迭代器到起始位置的距离，即索引
	int max_index = distance(cand_key_match.begin(), max_it);

	cout << "Maximum element: " << *max_it << endl;
	cout << "Index of maximum element: " << max_index << endl;

	return 0;
}