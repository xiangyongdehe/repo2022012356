## Part III. 3 ChIP-seq
## 致理-生21 孙铭一 2022012361
# 🧬 ChIP-seq Homework

## 1. 为什么在 ChIP-seq 实验中要做一个 control（通常叫 input）实验？

在 ChIP-seq 实验中，control（input）实验用于控制背景信号，其主要作用包括：

- **校正染色质开放性导致的非特异性富集**；
- **排除由建库、PCR 扩增、测序偏好等技术性噪音带来的误差**；
- **提高 peak calling 的特异性和准确性**；
- **提供与 IP 样本比较的基线**，从而识别真实的蛋白质-DNA 结合位点。

因此，input 是不可或缺的对照数据，确保结果的可靠性。

---

## 2. `findPeaks` 和 `findMotifsGenome.pl` 的主要参数解释

### 🟢 `findPeaks`
常用参数说明如下：
| 参数            | 含义                                      |
|-----------------|-------------------------------------------|
| `-style`        | 选择 peak 类型（`factor`：尖峰；`histone`：宽峰） |
| `-i`            | 输入 control 文件所在的 tagDirectory 目录 |
| `-F`            | 设置最小 Fold Change                      |
| `-P`            | 设置显著性阈值（p-value）                 |
| `-o auto`       | 自动生成输出文件名                        |

### 🟢 `findMotifsGenome.pl`
常用参数说明如下：
| 参数             | 含义                                 |
| -------------- | ---------------------------------- |
| `<peak file>`  | 输入 peak 文件（可以是 txt 或 bed）          |
| `<genome>`     | 基因组版本（如 `mm10`, `hg38`, `sacCer3`） |
| `<output dir>` | motif 输出文件夹                        |
| `-size`        | 设置 motif 提取区域大小，默认 200 bp          |
| `-len`         | 设置 motif 长度，支持多个（如 `8,10,12`）      |

---
## 3. 用homer重复本章中介绍的peak calling和motif finding分析
```bash
# 创建 Tag Directory（标签目录）
makeTagDirectory ~/chip-seq/output_homework/ip_dir ~/chip-seq/homework/ip.chrom_part.bam
makeTagDirectory ~/chip-seq/output_homework/input_dir ~/chip-seq/homework/input.chrom_part.bam

# 进行Peak Calling
findPeaks ~/chip-seq/output_homework/ip_dir \
  -style factor \
  -i ~/chip-seq/output_homework/input_dir \
  -o ~/chip-seq/output_homework/part.peak.txt \
  -F 8 -P 1e-8

# Motif分析
findMotifsGenome.pl ~/chip-seq/output_homework/part.peak.txt sacCer3 ~/chip-seq/output_homework/part.motif.output/ -size 200 -len 8,10,12
```
![image](https://github.com/user-attachments/assets/5fa3b963-a3e0-499d-ae29-4a5edb52c2ff)
![image](https://github.com/user-attachments/assets/3781d873-8838-4879-ae56-64f8a579cfaf)
