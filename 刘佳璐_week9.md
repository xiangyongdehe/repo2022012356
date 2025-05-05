#### 请解释在ChIP-seq实验中为什么一般都要平行做一个 control （通常叫 input）的实验。
##### 一、核心目的
**Input对照实验用于区分特异性结合信号与非特异性背景信号**，解决以下问题：
1. ​**排除开放染色质区域的非特异富集**  
   （DNA可及性高的区域更容易被随机剪切或沉淀）。
2. ​**校正DNA片段化偏好性**  
   （超声破碎时基因组的物理断裂存在区域偏好性）。
3. ​**消除实验技术性噪声**  
   （如测序偏差、抗体非特异性结合、PCR扩增偏好）。

##### 二、Input对照的具体作用
| ​**分析阶段**       | ​**功能说明**                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| ​**Peak Calling**    | 构建背景模型，识别显著富集区域（如MACS2使用λ参数建模）                     |
| ​**信号标准化**      | 提供基准水平，校正测序深度的差异（如TPM计算）                             |
| ​**假阳性过滤**      | 排除转录起始位点（TSS）、重复序列区等非特异高信号区域                     |
| ​**功能可视化**      | 在IGV等工具中直接比对IP与Input信号轨迹（如验证CTCF结合位点的特异性）      |


#### 请解释 findPeaks 和 findMotifsGenome.pl 主要参数的含义。
##### ​**1. `findPeaks` 核心参数**

| 参数         | 类型   | 作用                                   | 示例/推荐值             |
|--------------|--------|----------------------------------------|-------------------------|
| `-style`     | 必选   | 分析类型 (`factor`/`histone`)          | `-style factor`         |
| `-i`         | 必选   | Input对照的tag目录                     | `-i input_tagdir`       |
| `-o`         | 必选   | 输出文件名                             | `-o peaks.txt`          |
| `-F`         | 可选   | IP/Input的最小折叠变化阈值             | `-F 4` (默认=4)        |
| `-P`         | 可选   | p-value阈值 (泊松分布)                 | `-P 1e-5`              |
| `-size`      | 可选   | 转录因子的峰宽 (bp)                    | `-size 200` (默认=200) |
| `-minDist`   | 可选   | 相邻峰的最小间距                       | `-minDist 100`          |


##### ​**2. `findMotifsGenome.pl` 核心参数**

| 参数         | 类型   | 作用                                      | 示例/推荐值            |
|--------------|--------|-------------------------------------------|------------------------|
| peak文件     | 必选   | 输入的peak文件 (HOMER/BED格式)            | `peaks.txt`           |
| 基因组       | 必选   | 参考基因组版本                            | `sacCer3`/`hg38`      |
| 输出目录     | 必选   | 结果保存目录                              | `motif_output`        |
| `-len`       | 必选   | Motif长度 (单值或范围)                    | `-len 8,10`           |
| `-size`      | 必选   | 分析区域大小 (bp)                         | `-size 200`           |
| `-mask`      | 可选   | 屏蔽重复序列                              | `-mask`               |
| `-p`         | 可选   | 并行CPU核心数                             | `-p 4`                |
| `-bg`        | 可选   | 指定背景序列文件                          | `-bg input.txt`       |


#### 我们在容器的/home/test/chip-seq/homework目录中提供了酵母Snf1蛋白CHIP-seq的bam文件，ip.chrom_part.bam为IP实验数据，input.chrom\_part.bam为背景数据。请大家从这两个文件出发，用homer重复本章中介绍的peak calling和motif finding分析。请大家提交找到的motif的截图，以及Fold Change (vs Control) >=8且p-value (vs Control) < 10的−8次方的peaks(建议放在同一个文件中提交)。
##### step 1 运行 docker 并挂载容器
```bash
docker exec -it bioinfo_tsinghua bash
cd /home/test/chip-seq/
```


##### step 2 创建tag目录
```bash
makeTagDirectory ip_tagdir ip.chrom_part.bam
makeTagDirectory input_tagdir input.chrom_part.bam
```

##### step 3 peak calling
```bash
findPeaks ip_tagdir -style factor -i input_tagdir -o peaks0.txt
  -style factor \
  -i ~/chip-seq/output_homework/input_dir \
  -o ~/chip-seq/output_homework/peak0.txt \
  -F 8 -P 1e-8
```

##### step 4 motif分析
```bash
 findMotifsGenome.pl peaks0.txt sacCer2 motif_output0 -len 10
  ls -lh motif_output0 \
```
##### 输出结果
<img width="1109" alt="Screenshot 2025-05-05 at 20 46 36" src="https://github.com/user-attachments/assets/d1b8f254-f6be-4da4-a9b4-861e1bc2ba10" />


其中Fold Change (vs Control) >=8且p-value (vs Control) < 10的−8次方的peaks
<img width="1283" alt="Screenshot 2025-05-05 at 19 52 43" src="https://github.com/user-attachments/assets/4b8102d0-aa6d-4c83-8569-5dc98f40ca6a" />


