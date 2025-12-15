# 面向平方根无迹卡尔曼滤波器的硬件加速实现：集成Cholesky算子及其在PYNQ-Z2平台的系统设计

## Overview
本项目于2025年全国大学生嵌入式芯片与系统设计竞赛——AMD命题式赛道获得全国一等奖。
非常适合PYNQ的基础入门学习，代码量适中，结构清晰易调整
本项目基于 Xilinx HLS 和 Vivado，在 Zynq PYNQ-Z2 平台上实现平方根无迹卡尔曼滤波器（SR-UKF）的硬件加速。核心包括：
- 将 UKF 的时间更新和量测更新流程以硬件友好的方式重构
- 集成基于 Cholesky 分解的平方根形式协方差更新算子
- 通过 AXI4-Lite 接口将加速 IP 集成到 PS-PL 系统中，并在 PYNQ 上完成端到端系统设计与验证

## Project Structure

项目目录推荐组织形式如下：

```text
yolo_ukf_track/
├── notebooks/                   # 在 PYNQ‑Z2 上运行的 Jupyter Notebook
│   ├── yolo_ukf_soft.ipynb      # 纯软件版 YOLO + UKF 流程
│   ├── yolo_ukf_hardware.ipynb  # 调用 FPGA UKF IP 的硬件加速版
│   ├── design_1.bit             # 覆盖层比特流文件
│   ├── design_1.hwh             # 硬件描述文件（PYNQ 加载时使用）
│   ├── models/                  # YOLO 权重、配置及类别文件
│   ├── videos/                  # 测试视频（如 personwalking.avi 等）
│   └── result/                  # 检测/跟踪生成的中间结果和可视化视频
├── src/
│   ├── kernel/                  # HLS UKF 加速核工程
│   │   ├── ukf_accel.cpp        # HLS 顶层内核实现（ukf_accel_step 等）
│   │   ├── ukf_tb.cpp           # C 级测试平台，与 MATLAB 结果对比
│   │   ├── ukf.hpp              # SR‑UKF 核心算法与线性代数辅助函数
│   │   ├── cholesky.hpp         # Cholesky 分解/更新算子（XF Solver 改写）
│   │   ├── description.json     # HLS 内核描述文件
│   │   ├── hls_config.cfg       # HLS 构建配置
│   │   └── run_hls.tcl          # 一键运行 HLS 的 TCL 脚本
│   └── overlay/                 # Vivado 覆盖层工程
│       └── vivado_bd.tcl        # 创建 Block Design 以及Vivado运行的脚本
└── README.md                    # 项目说明与使用指南
```

## Hardware Requirements
- PYNQ-Z2 开发板（ XC7Z020 ）
- Micro-USB 电源及 USB 线缆
- 网络连接（用于访问 PYNQ Jupyter 环境）

## Software Requirements
- PYNQ 官方镜像（推荐与 PYNQ-Z2 匹配的稳定版本）
- Vivado 与 Vitis HLS 2024.2

## Getting Started

本节面向刚接触 PYNQ 与 HLS 的入门用户，给出从零开始的完整复现步骤。如在 Jupyter 上运行遇到“找不到模型/视频/结果目录”等问题，请先检查并修改 Notebook 中与路径相关的代码；视频的输入文件名需要根据你自己的数据进行设置。

Jupyter需要修改路径的代码位置：
- 纯软件版流程：`notebooks/yolo_ukf_soft.ipynb`
  - 路径与模型配置（根目录、模型目录等）：
    - `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:766-773` 中的 `SRC_DIR`、`MODELS_DIR`、`CFG_PATH`、`WEIGHTS_PATH`、`NAMES_PATH`、`DATA_PATH`
    - `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:1036-1038` 中的 `SRC_DIR`、`VIDEO_DIR`、`DEBUG_ROOT`
  - 输入视频所在目录：
    - `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:971` 与 `:1686` 中的 `video_dir = SRC_DIR`（默认从当前目录下的视频目录中批量读取）
- 硬件加速版流程：`notebooks/yolo_ukf_hardware.ipynb`
  - 路径与模型配置（根目录、模型目录等）：
    - `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_hardware.ipynb:174-181` 中的 `SRC_DIR`、`MODELS_DIR`、`CFG_PATH`、`WEIGHTS_PATH`、`NAMES_PATH`、`DATA_PATH`
  - 输入视频所在目录：
    - `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_hardware.ipynb:338` 与 `:805` 中的 `video_dir = SRC_DIR`
result_ukf_matlab.txt是Matlab源码进行批量测试数据的运行结果，Matlab代码在后面的`Performance`有提供地址

### 1. 在PC机上准备工程

1. 克隆或下载本项目到本地开发机：
   - 目录结构应包含 `yolo_ukf_track/notebook` 与 `yolo_ukf_track/src`。
2. 确认已安装软件环境：
   - Vivado 与 Vitis HLS 2024.2
   - 任意支持 `scp` 的终端工具（用于将文件拷贝到 PYNQ 板）

如果只想快速在 PYNQ 上体验硬件加速，而不修改硬件实现，可以直接使用仓库中已经生成好的：

- `notebook/design_1.bit`
- `notebook/design_1.hwh`

跳过“Building from Source”章节，直接看下面的“2. 部署到 PYNQ 板”。

### 2. 部署到 PYNQ 板并运行 Notebook

1. 启动 PYNQ‑Z2，确保网络连通，记下板子的 IP 地址（例如 `192.168.2.99`）。
2. 在开发机终端中，将整个 `yolo_ukf_track` 目录拷贝到板子上（示例命令）：

   ```bash
   scp -r yolo_ukf_track xilinx@<PYNQ_IP>:/home/xilinx/
   ```

   将 `<PYNQ_IP>` 替换为实际 IP，例如 `192.168.2.99`。

3. 在浏览器中访问：

   ```text
   http://<PYNQ_IP>:9090
   ```

   使用默认用户密码均是 `xilinx` 登录 Jupyter。

4. 在 Jupyter 文件浏览器中依次进入：

   ```text
   /home/xilinx/yolo_ukf_track/notebook
   ```
5. yolo部分已经进行抽帧检测，开发者可自行调整抽帧间隔以运行出更高帧率的视频
6. 运行纯软件版 YOLO+UKF 流程（不依赖 FPGA bitstream，仅做算法对齐）：
   - 打开 `yolo_ukf_soft.ipynb`
   - 从上到下依次执行各个 Cell
   - Notebook 会使用当前目录下的模型文件（`models/`）、示例视频（如 `personwalking.avi`），并在 `result/` 目录下生成检测与跟踪结果及输出视频

7. 运行硬件加速版 UKF（需要 bitstream）：
   - 打开 `yolo_ukf_hardware.ipynb`
   - 确认 `notebook/` 目录中存在：
     - `design_1.bit`
     - `design_1.hwh`
   - 依次运行 Notebook 单元：
     - Notebook 会在 PYNQ 板上加载 `design_1.bit`，通过 PYNQ 驱动 UKF 硬件 IP
     - YOLO 检测结果由 CPU 侧生成，UKF 状态更新在 PL 侧完成
     - 输出视频保存在 `notebook/result/<视频名>/out/` 目录（例如 `result/bird/out/video_result_yolov3.mp4`）

> 备注：软件版和硬件版的 UKF 算法是一一对应的，方便对比纯软件实现与硬件加速实现的轨迹和性能。

## Building from Source

本节说明如何从源码重新构建 HLS 内核与 Vivado 覆盖层。对初学者而言，可以完全按照以下命令逐步执行。

### HLS 构建步骤

1. 打开终端，进入 HLS 工程目录：

   ```bash
   cd yolo_ukf_track/src/kernel
   ```

2. 使用提供的 TCL 脚本运行 HLS 流程：

   ```bash
   vitis_hls run_hls.tcl
   ```

   默认配置中：
   - `CSIM` 置为 `1`，会执行 C 仿真，验证功能正确性
   - `CSYNTH`、`COSIM` 默认关闭（值为 `0`）

3. 如果希望继续执行综合与协同仿真，可以编辑 `run_hls.tcl` 顶部的开关，将对应变量改为 `1` 后重新运行：
   - `CSYNTH`：控制 C‑to‑RTL 综合
   - `COSIM`：控制 C/RTL 协同仿真

4. HLS 完成后，会在当前目录下生成工程：

   ```text
   ukf.prj/sol1/...
   ```

   其中包含：
   - HLS 报告（时序、资源占用）
   - 若启用了导出，包含生成的 IP（在后续 Vivado 覆盖层构建中使用）

5. 可以运行内置测试平台，检查 HLS 结果与参考结果的一致性：

   ```bash
   cd yolo_ukf_track/src/kernel
   g++ ukf_tb.cpp ukf_accel.cpp -I. -std=c++14 -o ukf_tb_host
   ./ukf_tb_host
   ```

   程序会读取 `inputs_ukf_common.txt`，生成 `result_ukf_hls.txt`，并与 `result_ukf_matlab.txt` 对比，输出误差指标与 PASS/FAIL 信息。

### Vivado 构建步骤

1. 打开终端，进入 overlay 目录：

   ```bash
   cd yolo_ukf_track/src/overlay
   ```

2. 使用提供的 TCL 脚本在批处理模式下运行 Vivado，自动创建 Block Design 并生成 bitstream：

   ```bash 
   vivado -mode batch -source vivado_bd.tcl
   ```

3. 脚本执行成功后，会在 `src/overlay/vivado_bd` 目录下生成：

   - `design_1.bit`
   - `design_1.hwh`注意：hwh需要解压一下，否则可能无法被识别
   - 以及中间工程文件

4. 将生成的文件拷贝到 Notebook 目录，供 PYNQ 加载使用（如需覆盖仓库中已有版本）：

   ```bash
   cp yolo_ukf_track/src/overlay/vivado_bd/design_1.bit  yolo_ukf_track/notebook/design_1.bit
   cp yolo_ukf_track/src/overlay/vivado_bd/design_1.hwh yolo_ukf_track/notebook/design_1.hwh
   ```

5. 重新按照 “Getting Started → 部署到 PYNQ 板并运行 Notebook” 的步骤，将更新后的 `yolo_ukf_track` 拷贝到 PYNQ 板，在 Jupyter 中重新运行 `yolo_ukf_hardware.ipynb`，即可在板端验证新生成的 bitstream。

## Performance
本项目在 PYNQ‑Z2 上的性能表现主要体现在：

- 使用硬件 UKF IP 后，单步滤波延时显著低于纯软件实现
- Notebook 中的示例会打印 UKF 调用耗时与平均耗时，可用于与软件版对比

更详细的时序和资源占用信息可以在：

- `yolo_ukf_track/src/kernel/ukf.prj/sol1/syn/report/` 中找到 HLS 报告
- `yolo_ukf_track/src/overlay/vivado_bd/` 中通过 Vivado 报告查看

## References
参考文献和资料
**论文出处**：[Unscented Kalman Filter 经典论文](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=a665183562768e29d87ce3073fbcde564ae00768)

**参考代码**：[sr-ukf MATLAB实现](https://github.com/JJHu1993/sr-ukf)

## License
本项目采用 Apache License 2.0 开源协议发布。

- 您可以自由地使用、修改和分发本项目代码（包括商业用途）。
- 在分发本项目或基于本项目的修改版本时，请保留原始版权声明和 Apache-2.0 协议文本。
- 本项目“按现状”提供，不对适用于特定目的、无缺陷等做任何保证，作者不对使用本项目造成的任何损失负责。

部分源码（例如 `src/kernel/cholesky.hpp` 及相关 XF Solver 代码）来自 Xilinx/AMD 官方示例，已在文件头部标明 Apache-2.0 版权与协议条款，请在再分发时一并保留这些说明。

Apache License 2.0 协议全文可参考：

- 英文原文：https://www.apache.org/licenses/LICENSE-2.0
