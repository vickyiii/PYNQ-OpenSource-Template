# PYNQ Open-Source Project Contribution Guide

## 目录 (Table of Contents)

- [项目简介 (Introduction)](#项目简介-introduction)
- [标准目录结构 (Standard Directory Structure)](#标准目录结构-standard-directory-structure)
- [文件要求 (File Requirements)](#文件要求-file-requirements)
- [提交流程 (Submission Process)](#提交流程-submission-process)
- [示例项目 (Example Projects)](#示例项目-example-projects)

---

## 项目简介 (Introduction)

本指南旨在帮助开发者以统一、可复现的方式贡献 PYNQ 开源项目。通过遵循标准化的目录结构和文件要求，我们可以：

- 📦 确保项目的可移植性和可复现性
- 🔧 简化项目部署和测试流程
- 📚 提高项目文档的一致性
- 🤝 促进社区协作和代码共享

This guide helps developers contribute PYNQ open-source projects in a unified, reproducible manner.

---

## 标准目录结构 (Standard Directory Structure)

所有提交的 PYNQ 项目必须遵循以下目录结构：

```
project_name/
├── notebook/
│   ├── *.ipynb              # Jupyter notebook 演示文件
│   ├── *.bit                # FPGA 比特流文件
│   ├── *.hwh                # 硬件描述文件
│   └── [其他资源文件]        # 权重文件、测试图片等
├── src/
│   ├── kernel_name_1/       # HLS 内核文件夹
│   │   ├── kernel.cpp       # 内核实现
│   │   ├── kernel_tb.cpp    # 测试平台
│   │   ├── kernel.hpp       # 头文件
│   │   ├── description.json # 内核描述
│   │   ├── hls_config.cfg   # HLS 配置
│   │   └── run_hls.tcl      # HLS 构建脚本
│   ├── kernel_name_2/       # 其他 HLS 内核
│   ├── .../
│   └── overlay/             # Vivado 覆盖层设计
│       ├── block_design_z2.tcl    # 块设计 TCL 脚本
│       └── run_vivado_z2.tcl      # Vivado 构建脚本
└── README.md                # 项目说明文档
```

---

## 文件要求 (File Requirements)

### 1. notebook/ 目录

**必需文件：**

- `*.ipynb` - 至少一个 Jupyter Notebook 演示文件
- `*.bit` - FPGA 比特流文件（与目标板匹配）（注意：参考项目未包含bitstream, 请务必提交时包含）
- `*.hwh` - 硬件描述文件

（注意：参考项目未包含bitstream 以及 硬件描述文件, 请务必提交时包含）

**可选文件：**

- 权重文件（如 `*.txt`, `*.npy`）
- 测试数据集
- 图片和可视化资源

**要求：**

- Notebook 应包含清晰的使用说明和示例, 中英文都可以，英文更好
- 比特流必须能在 PYNQ-Z2 板上直接运行
- 代码注释应使用英文（推荐）

### 2. src/ 目录

#### HLS 内核子目录

每个 HLS 内核必须包含以下文件：

| 文件名                   | 说明           | 必需 |
| ------------------------ | -------------- | ---- |
| `<kernel_name>.cpp`    | 内核实现源代码 | ✅   |
| `<kernel_name>_tb.cpp` | 测试平台代码   | ✅   |
| `<kernel_name>.hpp`    | 头文件         | ✅   |
| `description.json`     | 内核描述文件   | ✅   |
| `hls_config.cfg`       | HLS 配置文件   | ✅   |
| `run_hls.tcl`          | HLS 构建脚本   | ✅   |

**description.json 格式示例：**

```json
{
  "name": "kernel_name",
  "description": "Brief description of what this kernel does",
  "hls": {
    "config": "hls_config.cfg"
  }
}
```

**要求：**

- 代码应包含充分的注释
- 测试平台应能独立验证功能正确性
- HLS 配置应明确目标器件和时钟频率

#### overlay/ 子目录

必需文件：

- `block_design_z2.tcl` - Vivado Block Design设计脚本
- `run_vivado_z2.tcl` - Vivado 工程自动化构建脚本

**要求：**

- TCL 脚本应能完全自动化重建设计
- 明确指定目标 PYNQ 板型号（如 PYNQ-Z2 or KV260）
- 包含必要的 IP 配置和连接

### 3. README.md

**必需内容：**

```markdown
# 项目名称

## Overview
项目简介（中英文）

## Project Structure
目录结构说明

## Hardware Requirements
- 目标 PYNQ 板型号
- 其他硬件要求

## Software Requirements
- PYNQ 版本
- 其他依赖库

## Getting Started
### 快速开始步骤
1. 上传文件到 PYNQ 板
2. 打开 Jupyter Notebook
3. ...

## Building from Source
### HLS 构建步骤
### Vivado 构建步骤

## Performance
性能指标（如有）

## References
参考文献和资料

## License
开源协议
```

---

## 提交流程 (Submission Process)

## 方式一：仓库地址以及信息提交到表格中：

https://jsj.top/f/XIkmMx

### 方式二：通过分支提交到本仓库

### Step 1: 准备项目

1. **组织文件结构**

   ```bash
   # 创建标准目录结构
   mkdir -p your_project/notebook
   mkdir -p your_project/src/overlay
   ```
2. **整理 HLS 内核**

   ```bash
   # 为每个内核创建独立目录
   mkdir -p your_project/src/kernel_name
   ```
3. **验证完整性**

   - ✅ 所有必需文件已包含
   - ✅ 文件命名符合规范
   - ✅ 目录结构正确

### Step 2: 测试验证

1. **HLS 仿真测试**

   ```bash
   cd src/kernel_name
   vitis_hls -f run_hls.tcl
   # 验证 C Simulation 和 C/RTL Co-simulation 通过
   ```
2. **硬件测试**

   - 在 PYNQ 板上运行 Notebook
   - 验证功能正确性
   - 记录性能指标
3. **文档检查**

   - README.md 内容完整
   - 代码注释充分
   - 使用说明清晰

### Step 3: 提交项目

1. **Fork 仓库**

   ```bash
   # Fork PYNQ-OpenSource-Template 仓库
   ```
2. **创建项目分支**

   ```bash
   git checkout -b add-your-project-name
   ```
3. **添加项目文件**

   ```bash
   # 将项目复制到仓库根目录
   cp -r your_project/ PYNQ-OpenSource-Template/

   git add your_project/
   git commit -m "Add: your_project - brief description"
   git push origin add-your-project-name
   ```
4. **创建 Pull Request**

   - 提供项目简介
   - 说明主要功能和创新点
   - 附上测试结果截图

## 检查清单 (Checklist)

提交前请确认：

- [ ] 项目目录结构符合标准
- [ ] `notebook/` 包含 `.ipynb`, `.bit`, `.hwh` 文件
- [ ] 每个 HLS 内核包含完整的 6 个文件
- [ ] `overlay/` 包含 Vivado 构建脚本
- [ ] `README.md` 内容完整且格式正确
- [ ] HLS C Simulation 测试通过
- [ ] 在 PYNQ 板上实际测试通过
- [ ] 代码包含充分注释
- [ ] 提交信息清晰明确

---

## 常见问题 (FAQ)

### Q1: 如果我的项目只有一个 HLS 内核怎么办？

**A:** 仍然遵循相同结构，在 `src/` 下创建一个内核目录即可。

### Q2: 是否支持其他 PYNQ 板型号？

**A:** 是的！请在 overlay 脚本中明确指定板型号，一般是 KV260 以及 PYNQ-Z2，并在 README 中说明。

### Q3: 可以使用其他开源协议吗？

**A:** 可以，但请在 README.md 中明确声明协议类型。

### Q4: 如何处理大文件（如数据集）？

**A:** 建议使用下载链接，在 README 中提供获取方式，不要直接提交大文件。

---

**感谢您为 PYNQ 开源社区做出贡献！ Thank you for contributing to the PYNQ open-source community!** 🎉
