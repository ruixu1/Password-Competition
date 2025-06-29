# 项目名称
弱AES算法的设计与安全性分析

## 1. 项目简介
该项目要求设计一个最弱的8比特非线性S盒替换AES算法原有的S盒，保持其他部件不变，确实能够攻击的最高轮数。

## 2. 环境依赖
说明运行项目所需的语言版本、依赖库、工具链等：

- 操作系统：Linux
- 编程语言：C (支持 C11)，C++，Python 3.10+
- 外部依赖：
  - [Gurobi Optimizer](https://www.gurobi.com/)（用于 MILP）
  - Cadical（SAT 求解器）
  - OpenMPI（用于并行）

## 3. 各个文件功能说明
diff_word_search.cpp: 寻找r轮所有能达到最少活跃S盒个数的差分模式, 将得到的全部差分模式导出到"sol_word_all.txt"文件中
diff_sat.cpp: 导入"sol_word_all.txt"文件，并分别对每个差分模式作为一个限制条件生成一个SAT模型，然后调用Cadical计算可行解，如果有解，将仿射矩阵结果导出到"matrix_sol.txt", 将差分迹结果导出到"diff_sol.txt"(注: 最后一共有8组差分迹，它们导出文件分别是"diff_sol_{i}.txt", 0<=i<8)
diff_transform.py: 将"matrix_sol.txt"和"diff_sol_{i}.txt"中的内容转化为16进制的数组表示, 比如"matrix_sol.txt"每一行对应一个16进制元素，这样就可以直接在c文件中使用一个uint8_t matrix_sol[8]来表示这个矩阵
diff_key_recovery0.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_0.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_0.txt"中
diff_key_recovery1.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_1.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_1.txt"中
diff_key_recovery2.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_2.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_2.txt"中
diff_key_recovery3.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_3.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_3.txt"中
diff_key_recovery4.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_4.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_4.txt"中
diff_key_recovery5.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_5.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_5.txt"中
diff_key_recovery6.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_6.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_6.txt"中
diff_key_recovery7.c: 调用OpenMPI，使用传统差分分析方法使用"diff_sol_7.txt"的循环差分迹进行182轮的差分攻击，得到候选子密钥，导出到"cand_key_7.txt"中
diff_key_recovery_det.c: 调用OpenMPI，导入8个"cand_key_i.txt"，生成全部可能的候选子密钥，随机产生10个明文对验证得到唯一正确子密钥。
diff_key_recovery_check.c: 调用OpenMPI, 重复100次传统差分分析方法使用"diff_sol_0.txt"的循环差分迹进行102轮的差分攻击，看有多少次产生的候选子密钥集合中包括正确子密钥，检查我们选择的密钥对于差分分析成功概率的影响

## 4. 程序运行方法
diff_sat.cpp: gcc diff_word_search.c -o diff_word_search
./diff_word_search

diff_sat.cpp: g++ diff_sat.cpp -o diff_sat -I/mnt/e/cadical/src -L/mnt/e/cadical/build -lcadical -std=c++11(注：具体以自己Cadical按照目录为标准)
./diff_sat

diff_transform.py: python3 diff_transform.py

diff_key_recoveryi.c(0<=i<8): mpicc -O3 -march=native -fopenmp diff_key_recoveryi.c -o diff_key_recoveryi
mpicc -O3 -march=native -mavx2 -msse4.2 -ffast-math -funroll-loops -fopenmp diff_key_recoveryi.c -o diff_key_recoveryi
mpirun -np 100 --host localhost:100 ./diff_key_recoveryi

diff_key_recovery_det.c: mpicc -O3 -march=native -fopenmp diff_key_recovery_det.c -o diff_key_recovery_det
mpicc -O3 -march=native -mavx2 -msse4.2 -ffast-math -funroll-loops -fopenmp diff_key_recovery_det.c -o diff_key_recovery_det
mpirun -np 100 --host localhost:100 ./diff_key_recovery_det

diff_key_recovery_det.c: mpicc -O3 -march=native -fopenmp diff_key_recovery_check.c -o diff_key_recovery_check
mpicc -O3 -march=native -mavx2 -msse4.2 -ffast-math -funroll-loops -fopenmp diff_key_recovery_check.c -o diff_key_recovery_check
mpirun -np 100 --host localhost:100 ./diff_key_recovery_check
