# 1.请阐述 RNA-seq 中归一化基因表达值的几种基本计算方法。
# RNA-seq中基因表达值归一化的基本计算方法

在RNA-seq数据分析中，归一化（Normalization）用于消除技术偏差（如测序深度、基因长度等），使不同样本间的基因表达值具有可比性。以下是几种基本计算方法：

---

## 1. ​**RPKM/FPKM**
### 原理
- ​**RPKM** (Reads Per Kilobase per Million mapped reads) 和 ​**FPKM** (Fragments Per Kilobase per Million mapped reads) 是最早提出的归一化方法。
- 同时考虑 ​**测序深度**​（总reads数）和 ​**基因长度** 的影响。
- RPKM用于单端测序数据，FPKM用于双端测序数据。

### 公式
$$
\text{RPKM} = \frac{\text{基因的reads数}}{\left( \frac{\text{基因长度 (kb)} \times \text{总reads数}}{10^6} \right)}
$$
$$
\text{FPKM} = \frac{\text{基因的fragments数}}{\left( \frac{\text{基因长度 (kb)} \times \text{总fragments数}}{10^6} \right)}
$$

### 优缺点
- ​**优点**：简单直观，适用于样本内基因表达比较。
- ​**缺点**：不同样本的RPKM/FPKM总和可能不同，导致跨样本比较困难。

---

## 2. ​**TPM (Transcripts Per Million)**
### 原理
- 对RPKM/FPKM的改进，通过归一化使所有样本的TPM总和一致。
- 优先考虑基因长度，再调整测序深度。

### 公式
$$
\text{TPM} = \frac{\text{基因的reads数 / 基因长度 (kb)}}{\sum (\text{所有基因的reads数 / 基因长度})} \times 10^6
$$

### 优缺点
- ​**优点**：样本间TPM总和固定为百万，更便于跨样本比较。
- ​**缺点**：依然假设所有基因的表达量无全局差异，可能忽略样本间生物学差异。

---

## 3. ​**DESeq2的Median of Ratios**
### 原理
- DESeq2默认的归一化方法，假设大部分基因未发生差异表达。
- 通过计算基因的几何均值（geometric mean）和每个样本的Size Factor进行调整。

### 步骤
1. 计算每个基因在所有样本中的几何均值。
2. 对每个样本，计算基因表达值与几何均值的比值。
3. 取比值的中位数作为样本的Size Factor。
4. 原始reads数除以Size Factor得到归一化表达值。

### 公式
![image](https://github.com/user-attachments/assets/fa2edd2e-6a5c-4b08-9571-157cb0264dec)

### 优缺点
- ​**优点**：适用于差异表达分析，对低表达基因鲁棒。
- ​**缺点**：假设大部分基因无差异，当差异基因占比较高时可能失效。

---

## 4. ​**EdgeR的TMM (Trimmed Mean of M-values)**
### 原理
- 选择一组参考基因（非差异表达基因），计算样本间的缩放因子。
- 基于对数表达比（M值）和绝对表达量（A值）的修剪均值。

### 步骤
1. 选择一个样本作为参考。
2. 计算所有基因的M值（对数表达比）和A值（平均表达量）。
3. 去除M值极端和低表达的基因，计算剩余基因的加权均值作为缩放因子。

### 公式
$$
\text{TMM Factor} = \exp \left( \frac{\sum_{g \in G} w_g M_g}{\sum_{g \in G} w_g} \right)
$$
其中：
$$
w_g = \frac{\text{reads}_g}{\text{reads}_g + 1}
$$

### 优缺点
- ​**优点**：对样本间表达分布差异较大的数据更鲁棒。
- ​**缺点**：依赖参考基因的选择，需谨慎处理高差异表达数据。

---

## 5. ​**CPM (Counts Per Million)**
### 原理
- 仅考虑测序深度，不考虑基因长度。
- 适用于基因长度相近或无需长度校正的场景（如miRNA-seq）。

### 公式
$$
\text{CPM} = \frac{\text{基因的reads数}}{\text{总reads数}} \times 10^6
$$

### 优缺点
- ​**优点**：计算简单，适用于快速分析。
- ​**缺点**：忽略基因长度和组成偏差，不适合跨样本比较。

---

## 方法选择建议
| 方法              | 适用场景                          | 工具/软件       |
|-------------------|-----------------------------------|----------------|
| ​**TPM**           | 样本间基因表达水平比较            | StringTie, RSEM|
| ​**DESeq2/EdgeR**  | 差异表达分析（有生物学重复）      | DESeq2, EdgeR  |
| ​**RPKM/FPKM**     | 单样本内基因表达可视化            | Cufflinks      |
| ​**CPM**           | 长度相近的RNA（如miRNA）快速分析  | edgeR, limma   |


# 2.根据下述图片描述，填出对应选项: E B C

# 3.通过软件计算，判断给出文件shape02数据是来自哪一种sequencing protocols （strand nonspecific, strand specific - forward, strand specific - reverse)，并选择合适的参数计算shape02的read count matrix，给出AT1G09530基因(PIF3基因)上的counts数目。
```bash
#  进入 Docker 容器
docker exec -it bioinfo_featurecount /bin/bash
cd /home/test
/usr/local/bin/infer_experiment.py -r GTF/Arabidopsis_thaliana.TAIR10.34.bed -i bam/shape02.bam
```
#### 这两个分数接近相等(47.69% vs 49.16%)，说明数据是strand nonspecific​（非链特异性）测序

```bash
#计算 Read Count Matrix
/home/software/subread-2.0.3-source/bin/featureCounts \
  -T 4 \               # 使用 4 个线程
  -s 0 \               # 非链特异性（strand nonspecific）
  -p \                 # 针对 paired-end 数据
  -a GTF/Arabidopsis_thaliana.TAIR10.34.gtf \  # 注释文件
  -o /home/test/result/shape02_counts.txt \    # 输出路径
  bam/shape02.bam
# 提取 AT1G09530 (PIF3) 的 Counts
grep "AT1G09530" /home/test/result/shape02_counts.txt
```
#### 输出结果是
```bash
AT1G09530       1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1   3075768;3075768;3075768;3076401;3076401;3076401;3076459;3076459;3076459;3077173;3077173;3077173;3077173;3077378;3077378;3077378;3077378;3077378;3077378;3078346;3078346;3078346;3078346;3078346;3078346;3078545;3078545;3078545;3078545;3078545;3078545;3078843;3078843;3078843;3078843;3078843;3078843;3078984;3078984;3078984;3078984;3078984;3078984 3075852;3075852;3075852;3077286;3076808;3076748;3076808;3077286;3076748;3077286;3077286;3077286;3077286;3078257;3078257;3078257;3078257;3078257;3078257;3078453;3078453;3078453;3078453;3078453;3078453;3078610;3078610;3078610;3078610;3078610;3078610;3078908;3078908;3078908;3078908;3078908;3078908;3079544;3079544;3079544;3079654;3079654;3079654 +;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+;+   2762    86
```

#### counts数目为：2762

#  4. tumor-transcriptome-demo.tar.gz提供了结肠癌(COAD)，直肠癌(READ)和食道癌(ESCA)三种癌症各50个样本的bam文件用featureCount计算产生的结果。请大家编写脚本将这些文件中的counts合并到一个矩阵中(行为基因，列为样本), 计算logCPM的Z-score，并用 heatmap 展示，提供代码和heatmap。根据heatmap可视化的结果，你认为这三种癌症中哪两种癌症的转录组是最相似的?
```r
# 设置工作目录
setwd("~/Desktop/tumor-transcriptome-demo")

# 获取所有子目录中的.txt文件路径
files <- list.files(path = c("COAD", "ESCA", "READ"), 
                    pattern = "\\.txt$",
                    full.names = TRUE,
                    recursive = TRUE)

# 检查文件数量是否正确
if(length(files) != 150) warning(paste("找到", length(files), "个文件，预期150个"))

# 创建读取单个文件的函数
read_count_file <- function(file){
  # 读取数据，跳过注释行
  data <- read.table(file, header = TRUE, skip = 1)
  
  # 提取样本ID（从文件名）
  sample_id <- gsub(".txt$", "", basename(file))
  
  # 获取count列（通常是最后一列）
  count_col <- ncol(data)
  
  # 返回两列数据框
  data.frame(
    Geneid = data$Geneid,
    Count = data[, count_col],
    stringsAsFactors = FALSE
  )
}

# 并行读取所有文件（加快速度）
library(parallel)
count_list <- mclapply(files, read_count_file, mc.cores = detectCores())

# 给列表元素命名（使用样本ID）
names(count_list) <- gsub(".txt$", "", basename(files))

# 合并所有数据
merged_counts <- Reduce(
  function(x, y) merge(x, y, by = "Geneid", all = TRUE),
  count_list
)

# 转换为count矩阵
rownames(merged_counts) <- merged_counts$Geneid
count_matrix <- as.matrix(merged_counts[, -1])  # 移除Geneid列
colnames(count_matrix) <- names(count_list)

# 处理NA值（RNA-seq中通常将NA替换为0）
count_matrix[is.na(count_matrix)] <- 0

# 检查矩阵维度
dim(count_matrix)

library(edgeR)

# 创建DGEList对象
dge <- DGEList(counts = count_matrix)

# 计算CPM并取log2
logcpm <- cpm(dge, log = TRUE)

# 计算Z-score
z_score <- t(scale(t(logcpm)))

library(pheatmap)

# 创建样本分组信息
sample_groups <- data.frame(
  CancerType = factor(rep(c("COAD", "READ", "ESCA"), each = 50)),
  row.names = colnames(z_score)
)

# 绘制热图
pheatmap(z_score,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = sample_groups,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Z-score of logCPM values",
         color = colorRampPalette(c("blue", "white", "red"))(50))

# 保存PNG

png("heatmap.png", width=1200, height=1000, res=150)
pheatmap(z_score,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = sample_groups,
         main = "Z-score of logCPM values")
dev.off()

```
![heatmap](https://github.com/user-attachments/assets/7672713f-86f9-4b56-87a4-92eb7f0a1baa)

