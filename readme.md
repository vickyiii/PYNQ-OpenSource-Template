# Hardware Acceleration of Square Root Unscented Kalman Filter: Cholesky Operator Integration and System Design on PYNQ-Z2

## Overview

This project won the **First Prize** in the **2025 National University Embedded Chip and System Design Competition (AMD Proposition Track)**. It serves as an excellent introductory example for PYNQ learning, featuring moderate code volume, clear structure, and ease of modification.

This project implements hardware acceleration of the **Square Root Unscented Kalman Filter (SR-UKF)** on the **PYNQ-Z2** platform using **Vitis HLS** and **Vivado 2024.2**. Key features include:

- **Hardware-friendly reconstruction** of the UKF time update and measurement update processes.
- **Integration of Cholesky decomposition-based** square root covariance update operators.
- **End-to-end system design and verification** on PYNQ by integrating the acceleration IP into the PS-PL system via the **AXI4-Lite** interface.

## Project Structure

The recommended project directory structure is as follows:

```text
yolo_ukf_track/
├── notebooks/                    # Jupyter Notebooks to run on PYNQ-Z2
│   ├── yolo_ukf_soft.ipynb       # Pure software version of YOLO + UKF workflow
│   ├── yolo_ukf_hardware.ipynb   # Hardware accelerated version invoking FPGA UKF IP
│   ├── design_1.bit              # Overlay bitstream file
│   ├── design_1.hwh              # Hardware handoff file (used by PYNQ)
│   ├── models/                   # YOLO weights, configuration, and class files
│   ├── videos/                   # Test videos (e.g., personwalking.avi)
│   └── result/                   # Intermediate results and visualization videos
├── src/
│   ├── kernel/                   # HLS UKF acceleration kernel project
│   │   ├── ukf_accel.cpp         # HLS top-level kernel implementation (batch/single-step modes)
│   │   ├── ukf_tb.cpp            # C-level testbench, comparing with MATLAB results
│   │   ├── ukf.hpp               # SR-UKF core algorithms and linear algebra helpers
│   │   ├── cholesky.hpp          # Cholesky decomposition/update operator (adapted from XF Solver)
│   │   ├── description.json      # HLS kernel description file
│   │   ├── hls_config.cfg        # HLS build configuration
│   │   ├── inputs_ukf_common.txt # Input file for HLS testing
│   │   ├── result_ukf_matlab.txt # MATLAB test output results
│   │   └── run_hls.tcl           # TCL script for one-click HLS run
│   └── overlay/                  # Vivado overlay project
│       └── vivado_bd.tcl         # Script to create Block Design and run Vivado
└── README.md                     # Project documentation and user guide
```

## Hardware Requirements

- **PYNQ-Z2 Development Board** (XC7Z020)
- Micro-USB Power Supply and USB Cable
- Network Connection (for accessing PYNQ Jupyter environment)

## Software Requirements

- **PYNQ Official Image** (Recommended stable version matching PYNQ-Z2)
- **Vivado** and **Vitis HLS 2024.2**

## Getting Started

This section provides a complete reproduction guide from scratch for beginners new to PYNQ and HLS. If you encounter issues like "Model/Video/Result directory not found" when running on Jupyter, please check and modify the path-related code in the Notebook; the input video filename needs to be set according to your own data.

### Tips for Path and Data Preparation

Typical locations in Jupyter Notebooks that require path modification based on the local environment are as follows:

- **Pure Software Workflow**: `notebooks/yolo_ukf_soft.ipynb`
  - **Path and Model Configuration** (Root directory, model directory, etc.):
    - `SRC_DIR`, `MODELS_DIR`, `CFG_PATH`, `WEIGHTS_PATH`, `NAMES_PATH`, `DATA_PATH` in `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:766-773`
    - `SRC_DIR`, `VIDEO_DIR`, `DEBUG_ROOT` in `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:1036-1038`
  - **Input Video Directory**:
    - `video_dir = SRC_DIR` in `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_soft.ipynb:971` and `:1686` (Defaults to batch reading from the video directory in the current directory)
- **Hardware Acceleration Workflow**: `notebooks/yolo_ukf_hardware.ipynb`
  - **Path and Model Configuration** (Root directory, model directory, etc.):
    - `SRC_DIR`, `MODELS_DIR`, `CFG_PATH`, `WEIGHTS_PATH`, `NAMES_PATH`, `DATA_PATH` in `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_hardware.ipynb:174-181`
  - **Input Video Directory**:
    - `video_dir = SRC_DIR` in `/home/whp/Desktop/UKF/yolo_ukf_track/notebooks/yolo_ukf_hardware.ipynb:338` and `:805`

- `result_ukf_matlab.txt` contains the results of running batch test data with MATLAB source code. The link to the corresponding MATLAB reference implementation is provided in the **Performance** section later.

### 1. Prepare the Project on PC

1. **Clone or Download** this project to your local development machine:
   - The directory structure should include `yolo_ukf_track/notebooks` and `yolo_ukf_track/src`.
2. **Verify Installed Software Environment**:
   - Vivado and Vitis HLS 2024.2
   - Any terminal tool supporting `scp` (for copying files to the PYNQ board)

If you only want to quickly experience hardware acceleration on PYNQ without modifying the hardware implementation, you can directly use the pre-generated files in the repository:

- `notebook/design_1.bit`
- `notebook/design_1.hwh`

Skip the "Building from Source" section and go directly to "2. Deploy to PYNQ Board" below.

### 2. Deploy to PYNQ Board and Run Notebook

1. **Start PYNQ-Z2**, ensure network connectivity, and note the board's IP address (e.g., `192.168.2.99`).
2. **Copy the entire `yolo_ukf_track` directory** to the board using the terminal on your development machine (example command):

   ```bash
   scp -r yolo_ukf_track xilinx@<PYNQ_IP>:/home/xilinx/
   ```

   Replace `<PYNQ_IP>` with the actual IP, e.g., `192.168.2.99`.

3. **Access via Browser**:

   ```text
   http://<PYNQ_IP>:9090
   ```

   Login to Jupyter using the default username and password, both are `xilinx`.

4. **Navigate in Jupyter File Browser** to:

   ```text
   /home/xilinx/yolo_ukf_track/notebook
   ```

5. The YOLO part has already implemented frame skipping detection. Developers can adjust the frame skipping interval according to board performance and requirements to obtain a higher output frame rate while ensuring accuracy.

6. **Run Pure Software YOLO+UKF Workflow** (Does not depend on FPGA bitstream, only for algorithm alignment):
   - Open `yolo_ukf_soft.ipynb`
   - Execute each Cell sequentially from top to bottom
   - The Notebook will use model files in the current directory (`models/`), sample videos (e.g., `personwalking.avi`), and generate detection/tracking results and output videos in the `result/` directory

7. **Run Hardware Accelerated UKF** (Requires bitstream):
   - Open `yolo_ukf_hardware.ipynb`
   - Confirm the existence of:
     - `design_1.bit`
     - `design_1.hwh`
     in the `notebook/` directory.
   - Execute Notebook cells sequentially:
     - The Notebook will load `design_1.bit` on the PYNQ board and drive the UKF hardware IP via PYNQ.
     - YOLO detection results are generated on the CPU side, while UKF state updates are completed on the PL side.
     - Output videos are saved in the `notebook/result/<video_name>/out/` directory (e.g., `result/bird/out/video_result_yolov3.mp4`).

> **Note**: The software and hardware versions of the UKF algorithm correspond one-to-one, facilitating comparison of trajectories and performance between pure software implementation and hardware acceleration implementation.

## Building from Source

This section explains how to rebuild the HLS kernel and Vivado overlay from source code. Beginners can follow the commands below step by step.

### HLS Build Steps

1. **Open Terminal** and navigate to the HLS project directory:

   ```bash
   cd yolo_ukf_track/src/kernel
   ```

2. **Run HLS Flow** using the provided TCL script:

   ```bash
   vitis_hls run_hls.tcl
   ```

   In the default configuration:
   - `CSIM` is set to `1`, which executes C simulation to verify functional correctness.
   - `CSYNTH` and `COSIM` are disabled by default (values are `0`).

3. If you wish to proceed with synthesis and co-simulation, **edit the switches** at the top of `run_hls.tcl`, change the corresponding variables to `1`, and rerun:
   - `CSYNTH`: Controls C-to-RTL synthesis
   - `COSIM`: Controls C/RTL co-simulation

4. After HLS completion, a project will be generated in the current directory:

   ```text
   ukf.prj/sol1/...
   ```

   This includes:
   - HLS reports (timing, resource utilization)
   - Generated IP (if export is enabled), used in subsequent Vivado overlay construction

5. You can **run the built-in testbench** to check the consistency between HLS results and reference results:

   ```bash
   cd yolo_ukf_track/src/kernel
   g++ ukf_tb.cpp ukf_accel.cpp -I. -std=c++14 -o ukf_tb_host
   ./ukf_tb_host
   ```

   The program reads `inputs_ukf_common.txt`, generates `result_ukf_hls.txt`, compares it with `result_ukf_matlab.txt`, and outputs error metrics and PASS/FAIL information.

### Vivado Build Steps

1. **Open Terminal** and navigate to the overlay directory:

   ```bash
   cd yolo_ukf_track/src/overlay
   ```

2. **Run Vivado** in batch mode using the provided TCL script to automatically create the Block Design and generate the bitstream:

   ```bash 
   vivado -mode batch -source vivado_bd.tcl
   ```

3. Upon successful script execution, the following files will be generated in the `src/overlay/vivado_bd` directory:

   - `design_1.bit`
   - `design_1.hwh` (Note: If `design_1.hwh` is compressed, it must be decompressed first, otherwise PYNQ may not recognize it correctly)
   - Intermediate project files

4. **Copy the generated files** to the Notebook directory for PYNQ to load (if you need to overwrite the existing version in the repository):

   ```bash
   cp yolo_ukf_track/src/overlay/vivado_bd/design_1.bit  yolo_ukf_track/notebook/design_1.bit
   cp yolo_ukf_track/src/overlay/vivado_bd/design_1.hwh yolo_ukf_track/notebook/design_1.hwh
   ```

5. Follow the steps in "Getting Started → Deploy to PYNQ Board and Run Notebook" again to copy the updated `yolo_ukf_track` to the PYNQ board. Re-run `yolo_ukf_hardware.ipynb` in Jupyter to verify the newly generated bitstream on the board.

## Performance

The performance of this project on PYNQ-Z2 is highlighted by:

- **Significantly lower latency** for single-step filtering using the hardware UKF IP compared to pure software implementation.
- The Notebook example prints UKF execution time and average time, facilitating comparison with the software version.

More detailed timing and resource utilization information can be found in:

- **HLS Reports**: `yolo_ukf_track/src/kernel/ukf.prj/sol1/syn/report/`
- **Vivado Reports**: `yolo_ukf_track/src/overlay/vivado_bd/`

## References

**Paper**: [The Classical Unscented Kalman Filter Paper](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=a665183562768e29d87ce3073fbcde564ae00768)

**Reference Code**: [sr-ukf MATLAB Implementation](https://github.com/JJHu1993/sr-ukf)

## License

This project is licensed under the **Apache License 2.0**.

- You are free to use, modify, and distribute this project's code (including for commercial purposes).
- When distributing this project or modified versions based on it, please retain the original copyright notice and the Apache-2.0 license text.
- This project is provided "AS IS", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors be liable for any claim, damages or other liability.

Some source code (e.g., `src/kernel/cholesky.hpp` and related XF Solver code) is derived from official Xilinx/AMD examples and is marked with Apache-2.0 copyright and license terms in the file headers. Please retain these notices when redistributing.

The full text of the Apache License 2.0 can be found at:

- https://www.apache.org/licenses/LICENSE-2.0
