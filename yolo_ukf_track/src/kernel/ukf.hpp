#ifndef _XF_SOLVER_UKF_HPP_
#define _XF_SOLVER_UKF_HPP_

#include <cmath>
#include <cstring>
#include "hls_stream.h"
#include "ap_int.h"
#include "ap_fixed.h"
#include "cholesky.hpp"

// 本文件实现了一个基于平方根无迹卡尔曼滤波器（SR-UKF）的
// 小尺寸线性代数/滤波库，并对 Vitis HLS 做了适配。
// 约定：
// - N：状态维度（本工程中 N=3）
// - M：量测维度（本工程中 M=2）
// - x：状态均值向量
// - S：协方差矩阵的 Cholesky 因子（上三角），满足 P = S^T * S
// - q, r：过程噪声、量测噪声的标量超参数

namespace xf {
namespace solver {

// UKF 的权重参数结构体：
// - Wm：用于计算均值的权重
// - Wc：用于计算协方差的权重
// - cscale：协方差的缩放因子，用于生成 sigma 点
// - Wc_sqrt：对 Wc 取绝对值再开根号，用于平方根滤波中的协方差累积
template <int N>
struct UkfWeights {
    float Wm[2 * N + 1];
    float Wc[2 * N + 1];
    float cscale;
    float Wc_sqrt[2 * N + 1];
};

// 根据 Julier/UT 的公式计算权重。
// alpha、ki、beta 这些超参数控制 sigma 点分布的“散开”程度。
template <int N>
void ukf_compute_weights(UkfWeights<N>& w) {
    float alpha = 0.3f;
    float ki = 1.0f;
    float beta = 2.0f;
    float lambda = alpha * alpha * (N + ki) - N;
    float c = N + lambda;
    w.Wm[0] = lambda / c;
    for (int i = 1; i < 2 * N + 1; i++) w.Wm[i] = 0.5f / c;
    for (int i = 0; i < 2 * N + 1; i++) w.Wc[i] = w.Wm[i];
    w.Wc[0] = w.Wc[0] + (1 - alpha * alpha + beta);
    w.cscale = std::sqrt(c);
    for (int k = 0; k < 2 * N + 1; k++) w.Wc_sqrt[k] = std::sqrt(std::fabs(w.Wc[k]));
}

// 状态转移函数 f(x)：定义了系统的运动模型。
// 在本工程中，状态 x = [position, velocity, acceleration]^T 的一类非线性模型示例。
template <int N>
inline void f_state(const float x[N], float y[N]) {
    y[0] = x[1];
    y[1] = x[2];
    y[2] = 0.05f * x[0] * (x[1] + x[2]);
}

// 量测函数 h(x)：把状态映射到观测空间。
// 本例中是一个简单的“选前 M 个分量”的线性观测。
template <int N, int M>
inline void h_meas(const float x[N], float z[M]) {
    for (int i = 0; i < M; i++) z[i] = x[i];
}

// 在线性代数辅助函数部分，我们实现了一些小尺寸矩阵运算，
// 这些实现虽然不是通用 BLAS，但在 HLS 上更容易被综合与优化。

// 给矩阵 A 的对角线加上一个标量 d，用于增加过程噪声等。
template <int N>
inline void mat_add_diag(float A[N][N], float d) {
    for (int i = 0; i < N; i++) A[i][i] += d;
}

// P += w * v * v^T 的形式外积累加，用于协方差计算。
// 这里通过 UNROLL 指令展开循环，提高并行度。
template <int N, int K>
inline void outer_add_weighted(const float v[N], const float w, float P[N][N]) {
    for (int i = 0; i < N; i++) {
#pragma HLS UNROLL
        for (int j = 0; j < N; j++) {
#pragma HLS UNROLL
            P[i][j] += w * v[i] * v[j];
        }
    }
}

// 将矩阵 A 清零。
template <int N>
inline void mat_zero(float A[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) A[i][j] = 0.0f;
    }
}

// B = A 的逐元素拷贝。
template <int N>
inline void mat_copy(const float A[N][N], float B[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) B[i][j] = A[i][j];
    }
}

// C = A * B 的简单三重循环实现。
template <int N>
inline void mat_mul(const float A[N][N], const float B[N][N], float C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float s = 0.0f;
            for (int k = 0; k < N; k++) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
    }
}

// 支持指定有效 行/列 数的矩阵乘法，用于只有部分子矩阵有意义的场景。
template <int N>
inline void mat_mul_nm(const float A[N][N], const float B[N][N], int n, int m, float C[N][N]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            float s = 0.0f;
            for (int k = 0; k < n; k++) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
    }
}

// 单位矩阵 I 的构造。
template <int M>
inline void mat_identity(float I[M][M]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) I[j][i] = (i == j) ? 1.0f : 0.0f;
    }
}

// 使用高斯消元求逆矩阵（教学性质实现，尺寸较小）。
template <int M>
inline void mat_inverse_gauss(const float A[M][M], float invA[M][M]) {
    float aug[M][2 * M];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) aug[i][j] = A[i][j];
        for (int j = 0; j < M; j++) aug[i][M + j] = (i == j) ? 1.0f : 0.0f;
    }
    for (int col = 0; col < M; col++) {
        float piv = aug[col][col];
        float invp = 1.0f / piv;
        for (int j = 0; j < 2 * M; j++) aug[col][j] *= invp;
        for (int i = 0; i < M; i++) {
            if (i == col) continue;
            float factor = aug[i][col];
            for (int j = 0; j < 2 * M; j++) aug[i][j] -= factor * aug[col][j];
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) invA[i][j] = aug[i][M + j];
    }
}

// C = A * B，其中 A 为 N×M，B 为 M×M。
template <int N, int M>
inline void mat_mul_nm2(const float A[N][M], const float B[M][M], float C[N][M]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            float s = 0.0f;
            for (int k = 0; k < M; k++) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
    }
}

// 转置方阵 A -> AT。
template <int N>
inline void transpose(const float A[N][N], float AT[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) AT[j][i] = A[i][j];
    }
}

// 转置一般矩阵 A(N×M) -> AT(M×N)。
template <int N, int M>
inline void transpose_nm(const float A[N][M], float AT[M][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) AT[j][i] = A[i][j];
    }
}

// 生成 sigma 点：
// 给定当前均值 x 和协方差的 Cholesky 因子 S，以及权重 w，
// 生成 2N+1 个 sigma 点并存储在 X 中。
template <int N>
inline void make_sigma_points(const float x[N], const float S[N][N], const UkfWeights<N>& w,
                              float X[N][2 * N + 1]) {
    float A[N][N];
#pragma HLS PIPELINE II=1
    // A = cscale * S^T，用于构造正负偏移。
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            A[r][c] = w.cscale * S[c][r];
        }
    }
    // 第一个 sigma 点是均值本身。
    for (int i = 0; i < N; i++) X[i][0] = x[i];
#pragma HLS ARRAY_PARTITION variable=X complete dim=2
    // 后面 2N 个 sigma 点为 x +/- A 的列。
    for (int k = 0; k < N; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) {
            X[i][1 + k] = x[i] + A[i][k];
            X[i][1 + N + k] = x[i] - A[i][k];
        }
    }
}

// 无迹变换（UT）的“预测过程”部分：
// - 输入 sigma 点 X
// - 通过 f_state 传播，得到 Y
// - 对 Y 做加权得到预测均值 x1
// - 计算预测协方差的平方根 S1，并生成居中后的 sigma 点 X2
// - X1 保留了 f(x) 后未居中的 sigma 点，供后续量测 UT 使用
template <int N, int M>
inline void ukf_ut_process(const float X[N][2 * N + 1], const UkfWeights<N>& w, float q, float x1[N],
                           float S1[N][N], float X2[N][2 * N + 1], float X1[N][2 * N + 1]) {
    float Y[N][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=Y complete dim=2
#pragma HLS ARRAY_PARTITION variable=X2 complete dim=2
#pragma HLS ARRAY_PARTITION variable=X1 complete dim=2
#pragma HLS ARRAY_PARTITION variable=X complete dim=2
    // 1) 对每个 sigma 点做状态预测：Y = f(X)
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        float tmp[N];
        for (int i = 0; i < N; i++) tmp[i] = X[i][k];
        float ytmp[N];
        f_state<N>(tmp, ytmp);
        for (int i = 0; i < N; i++) Y[i][k] = ytmp[i];
    }
    // 2) 根据 Y 计算加权均值 x1。
    for (int i = 0; i < N; i++) x1[i] = 0.0f;
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) x1[i] += w.Wm[k] * Y[i][k];
    }
    // 3) X2 为居中后的 sigma 点，用于协方差计算。
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) X2[i][k] = Y[i][k] - x1[i];
    }
    // 4) X1 保留了未居中的 sigma 点（后续量测 UT 还会用到）。
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < N; i++) X1[i][k] = Y[i][k];
    }
    // 5) 使用平方根形式累积预测协方差：P = sum( w * v * v^T ) + q^2 I。
    float P[N][N];
    mat_zero<N>(P);
    for (int k = 0; k < 2 * N + 1; k++) {
        float wk = w.Wc_sqrt[k];
        float v[N];
        for (int i = 0; i < N; i++) v[i] = wk * X2[i][k];
        outer_add_weighted<N, 2 * N + 1>(v, 1.0f, P);
    }
    mat_add_diag<N>(P, q * q);
    // 6) 对 P 做 Cholesky 分解，得到新的平方根协方差 S1。
    float Lmat[N][N];
    xf::solver::choleskyTop<true, N, xf::solver::choleskyTraits<true, N, float, float>, float, float>(P, Lmat);
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) S1[r][c] = Lmat[r][c];
    }
}

// 无迹变换的“量测”部分：
// - 输入 X1：过程预测后的 sigma 点
// - 通过 h_meas 得到量测空间的 sigma 点 Z
// - 计算量测均值 z1 和量测协方差的平方根 S2
template <int N, int M>
inline void ukf_ut_meas(const float X1[N][2 * N + 1], const UkfWeights<N>& w, float r, float z1[M],
                        float S2[M][M], float Z2[M][2 * N + 1]) {
    float Z[M][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=Z complete dim=2
#pragma HLS ARRAY_PARTITION variable=Z2 complete dim=2
    // 1) 把状态空间的 sigma 点映射到量测空间。
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        float tmp[N];
        for (int i = 0; i < N; i++) tmp[i] = X1[i][k];
        float ztmp[M];
        h_meas<N, M>(tmp, ztmp);
        for (int i = 0; i < M; i++) Z[i][k] = ztmp[i];
    }
    // 2) 加权求量测均值。
    for (int i = 0; i < M; i++) z1[i] = 0.0f;
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < M; i++) z1[i] += w.Wm[k] * Z[i][k];
    }
    // 3) Z2 为居中后的量测 sigma 点。
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        for (int i = 0; i < M; i++) Z2[i][k] = Z[i][k] - z1[i];
    }
    // 4) 量测协方差累积，加上 r^2 噪声项。
    float P[M][M];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) P[i][j] = 0.0f;
    }
    for (int k = 0; k < 2 * N + 1; k++) {
        float wk = w.Wc_sqrt[k];
        float v[M];
        for (int i = 0; i < M; i++) v[i] = wk * Z2[i][k];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) P[i][j] += v[i] * v[j];
        }
    }
    for (int i = 0; i < M; i++) P[i][i] += r * r;
    // 5) Cholesky 分解得到量测协方差的平方根 S2。
    float Lmat[M][M];
    xf::solver::choleskyTop<true, M, xf::solver::choleskyTraits<true, M, float, float>, float, float>(P, Lmat);
    for (int r0 = 0; r0 < M; r0++) {
        for (int c0 = 0; c0 < M; c0++) S2[r0][c0] = Lmat[r0][c0];
    }
}

// 计算状态与量测之间的互协方差 P12。
template <int N, int M>
inline void cross_cov(const float X2[N][2 * N + 1], const float Z2[M][2 * N + 1], const UkfWeights<N>& w,
                      float P12[N][M]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) P12[i][j] = 0.0f;
    }
#pragma HLS ARRAY_PARTITION variable=X2 complete dim=2
#pragma HLS ARRAY_PARTITION variable=Z2 complete dim=2
    for (int k = 0; k < 2 * N + 1; k++) {
#pragma HLS PIPELINE II=1
        float wk = w.Wc[k];
        for (int i = 0; i < N; i++) {
#pragma HLS UNROLL
            for (int j = 0; j < M; j++) P12[i][j] += wk * X2[i][k] * Z2[j][k];
        }
    }
}

// 对上三角矩阵 R 做 rank-1 更新/降秩更新（cholupdate），
// 以保持平方根协方差的形式。downdate = true 表示减去 xx^T。
template <int N>
inline void cholupdate_upper(float R[N][N], const float x_in[N], bool downdate) {
    float x[N];
    for (int i = 0; i < N; i++) x[i] = x_in[i];
    float sign = downdate ? -1.0f : 1.0f;
    for (int k = 0; k < N; k++) {
        float rkk = R[k][k];
        float rkksq = rkk * rkk + sign * x[k] * x[k];
        if (rkksq < 0) rkksq = 0;
        float r = std::sqrt(rkksq);
        float c = r / rkk;
        float s = x[k] / rkk;
        R[k][k] = r;
        for (int j = k + 1; j < N; j++) {
            float rkj = R[k][j];
            float t = (rkj + sign * s * x[j]) / c;
            x[j] = c * x[j] - s * rkj;
            R[k][j] = t;
        }
    }
}

// UKF 的“更新”步骤：
// - 根据互协方差 P12 和量测协方差 S2 计算卡尔曼增益 K
// - 用观测 z 更新状态均值 x
// - 使用 cholupdate 对协方差的平方根进行 downdate，得到新的 S
template <int N, int M>
inline void ukf_update(const float x1[N], const float z[M], const float z1[M], const float S1[N][N], const float S2[M][M],
                       const float P12[N][M], float x[N], float S[N][N]) {
    float K[N][M];
    float U[N][M];
#pragma HLS ARRAY_PARTITION variable=K complete dim=2
#pragma HLS ARRAY_PARTITION variable=U complete dim=2
#pragma HLS ARRAY_PARTITION variable=S2 complete dim=2
#pragma HLS ARRAY_PARTITION variable=P12 complete dim=2
    float inv_diag[M];
#pragma HLS ARRAY_PARTITION variable=inv_diag complete dim=1
    // S2 为量测协方差的上三角矩阵，这里只用到了对角元素。
    for (int j = 0; j < M; j++) inv_diag[j] = 1.0f / S2[j][j];
    // 先做一次上三角前代，求解 U。
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        for (int j = 0; j < M; j++) {
            float s = 0.0f;
#pragma HLS UNROLL
            for (int k = 0; k < j; k++) s += S2[j][k] * U[i][k];
            U[i][j] = (P12[i][j] - s) * inv_diag[j];
        }
        // 再做一次回代，得到 K。
        for (int j = M - 1; j >= 0; j--) {
            float s2 = 0.0f;
#pragma HLS UNROLL
            for (int k = j + 1; k < M; k++) s2 += S2[k][j] * K[i][k];
            K[i][j] = (U[i][j] - s2) * inv_diag[j];
        }
    }
    // 创新项：实际量测减去预测量测。
    float innov[M];
    for (int i = 0; i < M; i++) innov[i] = z[i] - z1[i];
    // 更新状态均值。
    for (int i = 0; i < N; i++) {
        float s = 0.0f;
        for (int j = 0; j < M; j++) s += K[i][j] * innov[j];
        x[i] = x1[i] + s;
    }
    // 构造上三角矩阵 R = S1^T，后续对 R 做 cholupdate。
    float R[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) R[i][j] = S1[j][i];
    }
    // 对每一列量测方向做一次降秩更新，得到新的协方差平方根。
    for (int col = 0; col < M; col++) {
        float ucol[N];
        for (int i = 0; i < N; i++) ucol[i] = 0.0f;
        for (int k = 0; k < M; k++) {
            for (int i = 0; i < N; i++) ucol[i] += K[i][k] * S2[col][k];
        }
        cholupdate_upper<N>(R, ucol, true);
    }
    // 最终的平方根协方差 S 是 R 的转置。
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) S[i][j] = R[j][i];
    }
}

// 单步 SR-UKF：给定上一时刻状态（x_in, S_in）和本次量测 z，
// 输出更新后的状态（x_out, S_out）。这是硬件内核实际调用的核心函数。
template <int N, int M>
void ukf_step(const float z[M], float q, float r, float x_in[N], float S_in[N][N], float x_out[N], float S_out[N][N]) {
    UkfWeights<N> w;
    ukf_compute_weights<N>(w);
    float X[N][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=X complete dim=2
#pragma HLS DATAFLOW
    // 1) 生成 sigma 点。
    make_sigma_points<N>(x_in, S_in, w, X);
    float x1[N];
    float S1[N][N];
    float X2[N][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=X2 complete dim=2
    float X1p[N][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=X1p complete dim=2
    // 2) 过程无迹变换：传播到下一时刻，得到预测状态。
    ukf_ut_process<N, M>(X, w, q, x1, S1, X2, X1p);
    float z1[M];
    float S2[M][M];
    float Z2[M][2 * N + 1];
#pragma HLS ARRAY_PARTITION variable=Z2 complete dim=2
    // 3) 量测无迹变换：把 sigma 点映射到量测空间。
    ukf_ut_meas<N, M>(X1p, w, r, z1, S2, Z2);
    float P12[N][M];
    // 4) 计算互协方差并完成 Kalman 更新。
    cross_cov<N, M>(X2, Z2, w, P12);
    ukf_update<N, M>(x1, z, z1, S1, S2, P12, x_out, S_out);
}

// 多步运行接口：给定初始偏差 dx0 和一串噪声序列，
// 连续运行 N_steps 步 SR-UKF，返回最后一次的状态估计 x_out。
template <int N, int M, int MAX_STEPS>
void ukf_run(int N_steps, float q, float r, const float dx0[N], const float meas_noise[MAX_STEPS][M],
             const float proc_noise[MAX_STEPS][N], float x_out[N]) {
    // s：真实状态（仅在仿真环境中使用，在实际跟踪中通常来自系统动力学）
    float s[N];
    s[0] = 0.0f;
    s[1] = 0.0f;
    s[2] = 1.0f;
    // x：滤波器内部的状态估计，初始为真实状态加一个偏置 dx0。
    float x[N];
    for (int i = 0; i < N; i++) x[i] = s[i] + dx0[i];
    // S：协方差的平方根初始为单位矩阵。
    float S[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) S[i][j] = (i == j) ? 1.0f : 0.0f;
    }
    // 主循环：对每一步生成量测，调用 ukf_step 做一次 predict+update。
    for (int k = 0; k < N_steps; k++) {
        float zk[M];
        float hs[M];
        // 由真实状态生成理想量测，再叠加量测噪声。
        h_meas<N, M>(s, hs);
        for (int i = 0; i < M; i++) zk[i] = hs[i] + meas_noise[k][i];
        float x_next[N];
        float S_next[N][N];
        ukf_step<N, M>(zk, q, r, x, S, x_next, S_next);
        for (int i = 0; i < N; i++) x[i] = x_next[i];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) S[i][j] = S_next[i][j];
        }
        // 用状态方程推进真实状态，并叠加过程噪声。
        float fs[N];
        f_state<N>(s, fs);
        for (int i = 0; i < N; i++) s[i] = fs[i] + proc_noise[k][i];
    }
    // 返回最后一步的状态估计。
    for (int i = 0; i < N; i++) x_out[i] = x[i];
}

} // namespace solver
} // namespace xf

#endif

