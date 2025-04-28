# 1.多重检验校正与q值/p值的区别
## 什么是多重检验校正（Multiple Test Correction）？

### 背景
当同时进行大量统计检验（如基因差异表达分析中检测数万个基因）时，​**假阳性结果（False Positives）​**会显著增加。例如：
- 若对10,000个基因进行检验（α=0.05），即使所有基因均无差异，仍可能产生约500个假阳性结果（0.05×10,000）。

### 校正目的
通过调整显著性阈值（α）或p值，​**控制整体错误率**​（如假阳性比例）。

### 常用方法
| 方法                | 原理                                                                 | 特点                               |
|---------------------|----------------------------------------------------------------------|-----------------------------------|
| ​**Bonferroni**       | 调整α为α/m（m=检验次数）                                           | 严格，假阳性少但易漏真阳性（低功效） |
| ​**Benjamini-Hochberg (BH)** | 控制FDR（False Discovery Rate）                                     | 更灵活，适用于高通量数据           |
| ​**Holm**            | 逐步校正p值（比Bonferroni略宽松）                                   | 平衡严格性与功效                   |

---

## p值与q值的区别

### p值（P-value）
- ​**定义**：在零假设（H₀）成立时，观测到当前结果或更极端结果的概率。
- ​**解释**：  
  - 若p=0.01，表示在H₀成立时，有1%的概率出现此结果。
  - ​**控制错误类型**：Ⅰ类错误（假阳性）的概率。
- ​**局限**：未考虑多重检验问题，大量检验时假阳性累积。

### q值（q-value / FDR）
- ​**定义**：在拒绝零假设的显著结果中，​**假阳性比例（False Discovery Rate, FDR）​**的估计值。
- ​**解释**：  
  - 若q=0.05，表示在显著结果中，预期有5%是假阳性。
  - ​**控制错误类型**：错误发现率（FDR）。
- ​**应用场景**：高通量数据分析（如RNA-seq、GWAS）。

---

## 关键区别总结

| 特征                | p值                         | q值（FDR）                     |
|---------------------|-----------------------------|--------------------------------|
| ​**控制目标**         | 单次检验的Ⅰ类错误率         | 显著结果中的假阳性比例         |
| ​**调整方式**         | 未针对多重检验调整          | 基于多重检验校正（如BH方法）   |
| ​**解读示例**         | p=0.05 → 5%假阳性风险（单次）| q=0.05 → 显著结果中5%是假阳性  |
| ​**严格性**           | 更严格（单次检验）          | 更宽松（允许一定比例假阳性）   |
| ​**适用场景**         | 少量检验（如传统统计检验）   | 高通量数据（如基因组学）       |

---

## 示例说明
假设对10,000个基因进行差异表达分析：
- ​**未校正**：若设置p<0.05，可能有500个假阳性（5%×10,000）。
- ​**Bonferroni校正**：调整阈值为p<0.000005（0.05/10,000），假阳性≈0，但可能漏掉真实差异基因。
- ​**BH校正（FDR=0.05）​**：若筛选出500个显著基因，其中约25个是假阳性（5%×500），但能保留更多真实差异基因。

**优先选择**：  
- 探索性研究（需发现更多候选基因） → 使用q值（FDR）  
- 验证性研究（需严格避免假阳性） → 使用Bonferroni校正的p值

# DESeq2与edgeR的归一化方法详解

---

## 2.DESeq2的归一化方法：Median of Ratios

### ​**核心思想**
假设大部分基因在样本间**无差异表达**，通过计算样本间基因表达量的相对比值中位数，消除测序深度和RNA组成差异的影响。

---

### ​**具体步骤与公式**

#### 1. 计算基因的几何平均数
![image](https://github.com/user-attachments/assets/7e1587ce-58f4-42b1-b0ca-5b28d121b00a)


#### 2. 计算基因的比值
对每个样本i，基因g的比值是其计数与几何均值的比值：
<img width="363" alt="Screenshot 2025-04-28 at 21 38 42" src="https://github.com/user-attachments/assets/d49ad4d6-7c85-4dc9-95a6-d85cdfbb77aa" />


#### 3. 计算Size Factor
取样本i所有基因比值的中位数：
![image](https://github.com/user-attachments/assets/b9218016-cebe-4bbc-b23d-38f05c5c49ce)

#### 4. 归一化表达值
将原始read counts除以Size Factor：
![image](https://github.com/user-attachments/assets/2144b47f-6e69-42f9-a3c6-2c5931393e89)


---

### ​**关键特点**
- ​**依赖假设**：大部分基因无差异表达（适用于有生物学重复的数据）。  
- ​**鲁棒性**：中位数对异常值不敏感，适合低表达基因较多的数据。  
- ​**适用场景**：差异表达分析（需结合负二项分布模型）。

---

## 2.edgeR的归一化方法：TMM (Trimmed Mean of M-values)

### ​**核心思想**
基于**参考样本**​（通常选择表达量分布中等的样本），通过修剪极端值后的加权均值计算样本间的缩放因子，消除技术偏差。

---

### ​**具体步骤与公式**

#### 1. 选择参考样本
通常选择样本中位数或指定某个样本作为参考（记为样本 $r$）。

#### 2. 计算M值与A值
![image](https://github.com/user-attachments/assets/953093c5-6508-463c-b0b1-fb104358f16d)

#### 3. 修剪极端值
<img width="438" alt="Screenshot 2025-04-28 at 21 40 44" src="https://github.com/user-attachments/assets/5abc7fda-df1f-491d-810e-993b19a06978" />


#### 4. 计算TMM缩放因子
对保留的基因，计算加权平均M值，并转换为缩放因子：  
![image](https://github.com/user-attachments/assets/f8bf9eba-a2a0-41c1-b5cb-06b442df33fe)

#### 5. ​调整文库大小：
标准化后的文库大小为：
<img width="185" alt="Screenshot 2025-04-28 at 21 41 23" src="https://github.com/user-attachments/assets/452ba7c3-ff1b-4783-a9e4-09b376c551f7" />


### ​**关键特点**
- ​**灵活性**：通过修剪极端值提高鲁棒性，适合样本间RNA组成差异较大的数据。  
- ​**权重设计**：高表达基因的权重更大（信噪比更高）。  
- ​**适用场景**：差异表达分析（需结合经验贝叶斯方法）。

---

## 三、方法对比与选择建议

| 特征                | DESeq2 (Median of Ratios)       | edgeR (TMM)                     |
|---------------------|---------------------------------|---------------------------------|
| ​**核心假设**         | 大部分基因无差异表达            | 部分基因无差异表达（参考样本）  |
| ​**适用数据**         | 有生物学重复的样本              | 样本间RNA组成差异较大           |
| ​**鲁棒性**           | 对低表达基因敏感                | 对高表达基因更依赖（权重设计）   |
| ​**计算复杂度**       | 较高（需多次迭代）              | 较低（直接计算）                |
| ​**推荐场景**         | 常规差异分析（哺乳动物转录组）  | 特殊样本（如癌症vs正常组织）    |

---
# 3.利用我们以上介绍的方法和数据，分别使用DESeq2和edgeR找出uvr8突变型（uvr8）在光照前后的差异基因，保存为文本文件
```r
# 加载必要的包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "readr"))

library(DESeq2)
library(edgeR)
library(readr)

# 读取数据
count_data <- read.table("/Users/xiangyongdehe/Desktop/count_exon.txt", header = TRUE, row.names = 1)


# 检查数据维度
dim(count_data)

# 查看列名（样本名）
colnames(count_data)

# 查看样本数量
ncol(count_data)

# 检查数据维度
dim(count_data)

# 查看列名（样本名）
colnames(count_data)

# 查看样本数量
ncol(count_data)

# 提取样本信息
samples <- colnames(count_data)
genotype <- ifelse(grepl("^U", samples), "uvr8", "control")  # U开头是uvr8突变体
condition <- ifelse(grepl("1", samples), "light", "dark")     # 1代表light, 0代表dark

# 创建样本信息表
sample_info <- data.frame(
  genotype = factor(genotype),
  condition = factor(condition),
  row.names = samples
)

# 查看样本分组
print(sample_info)

#----DESeq2----
library(DESeq2)

# 只选择uvr8突变体的数据（UD开头的样本）
uvr8_samples <- grep("^UD", colnames(count_data), value = TRUE)
count_uvr8 <- count_data[, uvr8_samples]
sample_info_uvr8 <- sample_info[uvr8_samples, ]

# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_uvr8,
  colData = sample_info_uvr8,
  design = ~ condition
)

# 过滤低表达基因（至少3个样本的count > 10）
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 差异分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "light", "dark"))

# 保存结果
write.table(
  as.data.frame(res),
  file = "/Users/xiangyongdehe/Desktop/DESeq2_uvr8_light_vs_dark.txt",
  sep = "\t", quote = FALSE, col.names = NA
  
#----edge R 分析-----
library(edgeR)

# 只选择uvr8突变体的数据
y <- DGEList(counts = count_uvr8, group = sample_info_uvr8$condition)

# 过滤低表达基因
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM标准化
y <- calcNormFactors(y)

# 设计矩阵
design <- model.matrix(~ condition, data = sample_info_uvr8)

# 估计离散度
y <- estimateDisp(y, design)

# 差异分析
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)  # coef=2对应light vs dark

# 获取结果
edgeR_results <- topTags(qlf, n = Inf)

# 保存结果
write.table(
  edgeR_results$table,
  file = "/Users/xiangyongdehe/Desktop/edgeR_uvr8_light_vs_dark.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

```
### 结果txt文件在zip里

# 4.对于uvr8突变型的差异基因，定义|log2FC|>1，FDR<0.05的基因为差异表达基因。比较两个软件得到的差异基因有多少是重合的，有多少是不同的，用venn图的形式展示
```r
# 加载必要的包
library(DESeq2)
library(edgeR)
library(VennDiagram)
# 1. 数据准备（假设count_data和sample_info已加载）
uvr8_samples <- grep("^UD", colnames(count_data), value = TRUE)
count_uvr8 <- count_data[, uvr8_samples]
sample_info_uvr8 <- sample_info[uvr8_samples, ]
#  2. DESeq2分析
dds <- DESeqDataSetFromMatrix(countData = count_uvr8,
                              colData = sample_info_uvr8,
                              design = ~ condition)
dds <- DESeq(dds)
deseq2_res <- results(dds, contrast = c("condition", "light", "dark"))
# 3. edgeR分析
y <- DGEList(counts = count_uvr8, group = sample_info_uvr8$condition)
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data = sample_info_uvr8)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
edger_res <- topTags(qlf, n = Inf)$table
# 4. 提取显著差异基因
deseq2_sig <- rownames(deseq2_res)[abs(deseq2_res$log2FoldChange) > 1 & 
                                     deseq2_res$padj < 0.05 & 
                                     !is.na(deseq2_res$padj)]

edger_sig <- rownames(edger_res)[abs(edger_res$logFC) > 1 & 
                                   edger_res$FDR < 0.05 & 
                                   !is.na(edger_res$FDR)]
# 5. 统计结果
common_genes <- intersect(deseq2_sig, edger_sig)
deseq2_only <- setdiff(deseq2_sig, edger_sig)
edger_only <- setdiff(edger_sig, deseq2_sig)

# 6. 绘制Venn图
library(VennDiagram)

# 输入数据
deseq2_total <- 99
edger_total <- 50
overlap <- 36

# 计算特有基因数
deseq2_only <- deseq2_total - overlap  # 99-36=63
edger_only <- edger_total - overlap    # 50-36=14

# 绘制韦恩图
venn.plot <- draw.pairwise.venn(
  area1 = deseq2_total,    # DESeq2总基因数 (99)
  area2 = edger_total,     # edgeR总基因数 (50)
  cross.area = overlap,   # 重合基因数 (36)
  category = c("DESeq2", "edgeR"),
  fill = c("#1f78b4", "#e31a1c"),  # 蓝色和红色
  alpha = 0.6,
  lty = "blank",
  cex = 1.5,              # 数字大小
  cat.cex = 1.5,          # 类别标签大小
  cat.pos = c(-30, 30),   # 标签位置（角度）
  cat.dist = 0.05,        # 标签距离圆圈的距离
  margin = 0.1,
  print.mode = "raw"      # 显示实际数字
)

# 保存为高清PNG（600 dpi）
png("DEGs_Venn_uvr8.png", width = 2000, height = 2000, res = 300)
grid.draw(venn.plot)
dev.off()
```

![image](https://github.com/user-attachments/assets/c7cca2da-c427-46ea-ab57-b4702d121854)

# 5.对于edgeR找出的FDR<0.05的基因，选出log2FoldChange最大的10个基因和最小的10个基因。计算原始表达量矩阵的log10CPM值并对每个基因进行Z-score处理，使用刚才筛选出来的20个基因绘制热图（heatmap）作为最后结果输出。

```r
library(edgeR)
library(gplots)
library(RColorBrewer)

# 1. 数据准备 ----------------------------------------------------------------
# 假设已进行edgeR分析并得到edger_res结果
# 筛选FDR<0.05的显著基因
sig_genes <- edger_res[edger_res$FDR < 0.05 & !is.na(edger_res$FDR), ]

# 2. 选择logFC最大和最小的各10个基因 ----------------------------------------
top_10_up <- head(sig_genes[order(sig_genes$logFC, decreasing = TRUE), ], 10)
top_10_down <- head(sig_genes[order(sig_genes$logFC, decreasing = FALSE), ], 10)
selected_genes <- rbind(top_10_up, top_10_down)

# 3. 计算log10CPM并进行Z-score标准化 ---------------------------------------
# 假设y是已创建的DGEList对象（包含标准化因子）
logcpm <- cpm(y, log = TRUE, prior.count = 1)  # 计算log10CPM
logcpm_selected <- logcpm[rownames(selected_genes), ]

# Z-score标准化（按行）
zscore_matrix <- t(scale(t(logcpm_selected)))

# 4. 准备热图注释 ----------------------------------------------------------
# 假设sample_info_uvr8包含condition列（light/dark）
condition_colors <- ifelse(sample_info_uvr8$condition == "light", "gold", "darkblue")

# 5. 绘制红蓝经典热图 ------------------------------------------------------
pdf("edgeR_top20_genes_heatmap.pdf", width = 8, height = 10)

heatmap.2(
  zscore_matrix,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),  # 红蓝渐变色
  scale = "none",          # 已做Z-score无需再缩放
  trace = "none",         # 关闭trace线
  margins = c(10, 12),     # 下边距10，左边距12
  ColSideColors = condition_colors,  # 样本分组颜色条
  dendrogram = "row",     # 仅对基因聚类
  Rowv = TRUE,            # 对行聚类
  Colv = FALSE,           # 不对列聚类
  cexRow = 0.8 + 1/log10(nrow(zscore_matrix)),  # 自适应行标签大小
  cexCol = 0.8 + 1/log10(ncol(zscore_matrix)),  # 自适应列标签大小
  key = TRUE,             # 显示颜色标尺
  keysize = 1,            # 标尺大小
  density.info = "none",  # 关闭密度图
  main = "Top 20 DEGs in uvr8 (log10CPM Z-scores)",  # 标题
  labCol = sub("_", "\n", colnames(zscore_matrix))  # 处理样本名中的下划线
)

# 添加图例
legend("topright",
       legend = c("Light", "Dark"),
       fill = c("gold", "darkblue"),
       border = NA,
       bty = "n",
       cex = 0.8)

dev.off()

_____

```
## 热图结果

<img width="793" alt="Screenshot 2025-04-28 at 22 53 21" src="https://github.com/user-attachments/assets/cdc99989-8c69-4b7a-8711-c9c30e224384" />
