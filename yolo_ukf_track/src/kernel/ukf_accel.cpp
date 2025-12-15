#include "ukf.hpp"

// 该文件给出 HLS 工程的顶层 C 接口函数。
// 在综合时，Vitis HLS 会把这里的函数作为 PL 端 IP 的入口。

extern "C" {

// 多步运行版本：在 PL 中连续执行 N_steps 次 UKF 预测+更新，
// 输入为噪声序列，输出为最后一步的状态估计 x_out。
void ukf_accel(
    int N_steps,
    float q,
    float r,
    const float dx0[3],
    const float meas_noise[256][2],
    const float proc_noise[256][3],
    float x_out[3]) {
#pragma HLS INLINE off // 保持该函数为独立顶层，不内联到内部实现
    xf::solver::ukf_run<3, 2, 256>(N_steps, q, r, dx0, meas_noise, proc_noise, x_out);
}

// 单步版本：硬件中常用的接口形式。
// - 输入 z：当前时刻的量测（例如来自摄像头检测结果）
// - 输入 q, r：过程噪声、量测噪声的标量超参数
// - 输入 x_in, S_in：上一次的状态均值和协方差的 Cholesky 因子
// - 输出 x_out, S_out：本次更新后的状态和协方差的 Cholesky 因子
void ukf_accel_step(
    const float z[2],
    float q,
    float r,
    float x_in[3],
    float S_in[3][3],
    float x_out[3],
    float S_out[3][3]) {
#pragma HLS INLINE off
    // AXI4-Lite 控制接口声明：
    // 所有标量/小数组都通过控制总线读写，方便在 PYNQ 上由 CPU 配置。
#pragma HLS INTERFACE s_axilite port=z bundle=control
#pragma HLS INTERFACE s_axilite port=q bundle=control
#pragma HLS INTERFACE s_axilite port=r bundle=control
#pragma HLS INTERFACE s_axilite port=x_in bundle=control
#pragma HLS INTERFACE s_axilite port=S_in bundle=control
#pragma HLS INTERFACE s_axilite port=x_out bundle=control
#pragma HLS INTERFACE s_axilite port=S_out bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // 调用 UKF 的核心单步更新函数，模板参数 <3,2> 表示：
    // - 状态维度 N = 3（位置、速度等）
    // - 量测维度 M = 2（例如 x、y 位置）
    xf::solver::ukf_step<3, 2>(z, q, r, x_in, S_in, x_out, S_out);
}

}

