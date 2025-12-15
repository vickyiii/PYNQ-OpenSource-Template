#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "ukf.hpp"

// 这是 UKF 硬件内核的 C 级仿真/对比测试程序。
// 主要功能：
// 1. 从文本文件中读取统一的测试输入（初始偏差、过程噪声、量测噪声）
// 2. 使用 C 参考实现 xf::solver::ukf_step 生成“参考结果” result_ukf_hls.txt
// 3. 从 MATLAB 版本导出的 result_ukf_matlab.txt 中读取结果
// 4. 计算两者之间的误差指标（RMSE / MAE / NRMSE 等），检查是否通过阈值

extern "C" void ukf_accel(int N_steps, float q, float r, const float dx0[3], const float meas_noise[256][2],
                           const float proc_noise[256][3], float x_out[3]);
extern "C" void ukf_accel_step(const float z[2], float q, float r, float x_in[3], float S_in[3][3], float x_out[3],
                                float S_out[3][3]);

int main() {
    // 1. 打开测试输入文件，里面包含：
    //    - n, m, N：状态维度、量测维度、步数
    //    - q, r：过程噪声、量测噪声
    //    - 初始偏差 dx0
    //    - 每一步的量测噪声 meas_noise 与过程噪声 proc_noise
    const char* inpath = "inputs_ukf_common.txt";
    FILE* fp = fopen(inpath, "r");
    if (!fp) {
        // 兼容在 HLS 工程目录下运行的相对路径。
        inpath = "../../../../inputs_ukf_common.txt";
        fp = fopen(inpath, "r");
    }
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", inpath);
        return -1;
    }

    int n, m, N;
    if (fscanf(fp, "%d,%d,%d\n", &n, &m, &N) != 3) {
        fprintf(stderr, "ERROR: bad header\n");
        return -1;
    }

    float q, r;
    if (fscanf(fp, "%f,%f\n", &q, &r) != 2) {
        fprintf(stderr, "ERROR: bad noise line\n");
        return -1;
    }

    // 读取初始偏差 dx0（长度为 n）。
    float dx0[3] = {0, 0, 0};
    for (int i = 0; i < n; i++) {
        float v;
        if (fscanf(fp, i < n - 1 ? "%f," : "%f\n", &v) != 1) {
            fprintf(stderr, "ERROR: bad dx0 line\n");
            return -1;
        }
        dx0[i] = v;
    }

    // 每一步的量测噪声和过程噪声。根据输入文件的约定，最多 256 步。
    float meas_noise[256][2];
    float proc_noise[256][3];
    for (int k = 0; k < N; k++) {
        float z0, z1, pn0, pn1, pn2;
        if (fscanf(fp, "%f,%f,%f,%f,%f\n", &z0, &z1, &pn0, &pn1, &pn2) != 5) {
            fprintf(stderr, "ERROR: bad step line\n");
            return -1;
        }
        meas_noise[k][0] = z0;
        meas_noise[k][1] = z1;
        proc_noise[k][0] = pn0;
        proc_noise[k][1] = pn1;
        proc_noise[k][2] = pn2;
    }
    fclose(fp);

    // 2. 初始化真实状态 s 和滤波器内部状态 (x, S)。
    //    s 用于生成“理想量测”，x 为估计值，S 为协方差的平方根。
    float s[3] = {0.0f, 0.0f, 1.0f};
    float x[3] = {s[0] + dx0[0], s[1] + dx0[1], s[2] + dx0[2]};
    float S[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) S[i][j] = (i == j) ? 1.0f : 0.0f;
    }

    // 3. 打开输出文件，用于保存 C 参考实现的结果（与 MATLAB 结果对比）。
    FILE* fpout = fopen("result_ukf_hls.txt", "w");
    if (!fpout) {
        fprintf(stderr, "ERROR: cannot open result file for write\n");
        return -1;
    }
    fprintf(fpout, "%% Step, Actual[0], Actual[1], Actual[2], Estimate[0], Estimate[1], Estimate[2]\n");

    // 4. 主仿真循环：使用 xf::solver::ukf_step 进行 N 步滤波。
    for (int k = 0; k < N; k++) {
        // 由真实状态 s 生成理想量测，再叠加量测噪声，得到实际观测 zk。
        float zk[2];
        float hs[2];
        xf::solver::h_meas<3,2>(s, hs);
        zk[0] = hs[0] + meas_noise[k][0];
        zk[1] = hs[1] + meas_noise[k][1];

        // 调用 C 版本的 UKF 单步更新。
        float x_next[3];
        float S_next[3][3];
        xf::solver::ukf_step<3,2>(zk, q, r, x, S, x_next, S_next);

        // 把当前步的真实状态和估计结果写入文件，用于后处理/画图。
        fprintf(fpout, "%d, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
                k+1, s[0], s[1], s[2], x_next[0], x_next[1], x_next[2]);

        // 更新滤波器内部状态。
        for (int i = 0; i < 3; i++) x[i] = x_next[i];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) S[i][j] = S_next[i][j];
        }

        // 用状态方程推进真实状态，并加入过程噪声。
        float fs[3];
        xf::solver::f_state<3>(s, fs);
        for (int i = 0; i < 3; i++) s[i] = fs[i] + proc_noise[k][i];
    }

    // 5. 用最后一步的状态再运行一次硬件接口版本 ukf_accel_step，
    //    目的主要是检查硬件顶层接口和 C 版本的一致性。
    float hs_last[2];
    xf::solver::h_meas<3,2>(s, hs_last);
    float z_last[2];
    z_last[0] = hs_last[0] + meas_noise[N-1][0];
    z_last[1] = hs_last[1] + meas_noise[N-1][1];
    float x_out_step[3];
    float S_out_step[3][3];
    ukf_accel_step(z_last, q, r, x, S, x_out_step, S_out_step);
    fclose(fpout);

    // 6. 读取 MATLAB 结果文件和当前 C 结果文件，进行数值误差对比。
    const char* mpath = "result_ukf_matlab.txt";
    const char* hpath = "result_ukf_hls.txt";
    FILE* fm = fopen(mpath, "r");
    if (!fm) {
        // 同样兼容在 HLS 工程子目录下运行的相对路径。
        mpath = "../../../../result_ukf_matlab.txt";
        fm = fopen(mpath, "r");
    }
    FILE* fh = fopen(hpath, "r");
    if (!fm || !fh) {
        fprintf(stderr, "ERRORCHECK: cannot open results\n");
        if (fm) fclose(fm);
        if (fh) fclose(fh);
        return -2;
    }

    // 按行解析两个结果文件，分别存储步号和对应的估计值。
    char buf[512];
    int m_steps[4096];
    float m_est[4096][3];
    int h_steps[4096];
    float h_est[4096][3];
    int mc = 0, hc = 0;
    while (fgets(buf, sizeof(buf), fm)) {
        if (buf[0] == '%') continue;
        int st; float a0,a1,a2,e0,e1,e2;
        if (sscanf(buf, "%d, %f, %f, %f, %f, %f, %f", &st, &a0, &a1, &a2, &e0, &e1, &e2) == 7) {
            m_steps[mc] = st;
            m_est[mc][0] = e0; m_est[mc][1] = e1; m_est[mc][2] = e2;
            mc++;
        }
    }
    while (fgets(buf, sizeof(buf), fh)) {
        if (buf[0] == '%') continue;
        int st; float a0,a1,a2,e0,e1,e2;
        if (sscanf(buf, "%d, %f, %f, %f, %f, %f, %f", &st, &a0, &a1, &a2, &e0, &e1, &e2) == 7) {
            h_steps[hc] = st;
            h_est[hc][0] = e0; h_est[hc][1] = e1; h_est[hc][2] = e2;
            hc++;
        }
    }
    fclose(fm);
    fclose(fh);

    // 7. 逐步对齐步号，计算三维状态上的误差指标。
    double sum_err_sq[3] = {0.0,0.0,0.0};
    double sum_mat_sq[3] = {0.0,0.0,0.0};
    double sum_err_abs[3] = {0.0,0.0,0.0};
    double max_err[3] = {0.0,0.0,0.0};
    int cnt = 0;
    for (int i = 0; i < mc; i++) {
        int st = m_steps[i];
        for (int j = 0; j < hc; j++) {
            if (h_steps[j] == st) {
                for (int d = 0; d < 3; d++) {
                    double e = (double)h_est[j][d] - (double)m_est[i][d];
                    sum_err_sq[d] += e*e;
                    sum_err_abs[d] += std::fabs(e);
                    if (std::fabs(e) > max_err[d]) max_err[d] = std::fabs(e);
                    double mm = (double)m_est[i][d];
                    sum_mat_sq[d] += mm*mm;
                }
                cnt++;
                break;
            }
        }
    }
    if (cnt == 0) {
        fprintf(stderr, "ERRORCHECK: no common steps\n");
        return -3;
    }

    // 8. 汇总得到 RMSE / MAE / NRMSE 和聚合指标，用于 PASS/FAIL 判定。
    double rmse[3], nrmse[3], mae[3];
    for (int d = 0; d < 3; d++) {
        rmse[d] = std::sqrt(sum_err_sq[d] / (double)cnt);
        double denom = std::sqrt((sum_mat_sq[d] / (double)cnt) + 1e-12);
        nrmse[d] = rmse[d] / denom;
        mae[d] = sum_err_abs[d] / (double)cnt;
    }
    double agg_rmse = std::sqrt((sum_err_sq[0]+sum_err_sq[1]+sum_err_sq[2]) / (3.0*(double)cnt));
    double agg_denom = std::sqrt(((sum_mat_sq[0]+sum_mat_sq[1]+sum_mat_sq[2]) / (3.0*(double)cnt)) + 1e-12);
    double agg_nrmse = agg_rmse / agg_denom;
    printf("ERRORCHECK: cnt=%d\n", cnt);
    printf("ERRORCHECK: RMSE: %.8f, %.8f, %.8f\n", rmse[0], rmse[1], rmse[2]);
    printf("ERRORCHECK: MAE: %.8f, %.8f, %.8f\n", mae[0], mae[1], mae[2]);
    printf("ERRORCHECK: MaxAbs: %.8f, %.8f, %.8f\n", max_err[0], max_err[1], max_err[2]);
    printf("ERRORCHECK: NRMSE: %.8f, %.8f, %.8f\n", nrmse[0], nrmse[1], nrmse[2]);
    printf("ERRORCHECK: AggNRMSE: %.8f\n", agg_nrmse);

    // 设置通过阈值：如果归一化 RMSE 过大，则认为实现有问题。
    double thr = 0.30;
    if (agg_nrmse > thr) {
        fprintf(stderr, "ERRORCHECK: FAIL (AggNRMSE=%.6f > %.6f)\n", agg_nrmse, thr);
        return -4;
    }
    printf("ERRORCHECK: PASS\n");
    return 0;
}
