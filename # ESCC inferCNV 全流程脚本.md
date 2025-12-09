# ESCC inferCNV 全流程脚本（严格复现文章 2.4 方法参数）

✨ Overview

本脚本严格遵循 ESCC 单细胞文章第 2.3–2.4 部分，包括：
	•	Seurat → 上皮细胞提取
	•	20% EN 上皮作为 baseline
	•	cutoff = 0.1
	•	noise filter = 0
	•	基因按染色体排序
	•	每条染色体 101 基因滑动窗口平滑
	•	expression truncation [-3, 3]
	•	CNV score = CNV 偏移值平方和
	•	最终输出 CNV score，可用于 malignant/non-malignant 分类

⚠ 运行前，请手动修改路径、meta 列名、gene_order 文件路径等。

# 1️⃣ 数据准备（从 Seurat 出发）
```r
library(Seurat)
library(dplyr)
library(infercnv)
library(zoo)
set.seed(123)

## ---- 修改你的列名 ----
sample_col   <- "sample_type"   # ET / EN
celltype_col <- "celltype"      # Epithelial / Immune / Fibroblast...

meta <- seu@meta.data

# 只保留上皮细胞
meta_epi <- meta[meta[[celltype_col]] == "Epithelial", ]

epi_EN_cells <- rownames(meta_epi[meta_epi[[sample_col]] == "EN", ])
epi_ET_cells <- rownames(meta_epi[meta_epi[[sample_col]] == "ET", ])

## ---- 随机抽取 20% 作为 baseline ----
n_ref <- ceiling(0.2 * length(epi_EN_cells))
ref_cells <- sample(epi_EN_cells, size = n_ref, replace = FALSE)

query_cells <- epi_ET_cells
cells_use <- c(ref_cells, query_cells)
```

# 2️⃣ 输出 inferCNV 所需输入文件(这一步应该也不用做，先理解一下）
2.1 表达矩阵（好像是需要从seurat提取的）
```r
expr <- as.matrix(seu@assays$RNA@counts[, cells_use])

dir.create("infercnv_input", showWarnings = FALSE)
counts_file <- "infercnv_input/ESCC_epi_counts.txt"

write.table(expr, counts_file,
            sep = "\t", quote = FALSE, col.names = NA)
```
2.2 注释文件
```r
annot <- data.frame(
  cell  = cells_use,
  group = ifelse(cells_use %in% ref_cells, "Epi_EN_ref", "Epi_ET_query")
)

annot_file <- "infercnv_input/ESCC_epi_annotations.txt"
write.table(annot, annot_file,
            sep = "\t", quote = FALSE, row.names = FALSE)
```
2.3  基因位置信息（gene_order_file）
格式需包含四列： gene, chr, start, end必须与表达矩阵行名一致。
以下为示例路径
```r
gene_order_file <- "infercnv_input/gene_order_hg38.txt"
```
# 【重点！！】3️⃣ 创建 inferCNV 对象 + 运行（cutoff=0.1）
```r
infer_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_file,
  annotations_file  = annot_file,
  delim             = "\t",
  gene_order_file   = gene_order_file,
  ref_group_names   = c("Epi_EN_ref")
)

out_dir <- "infercnv_ESCC_epi"
dir.create(out_dir, showWarnings = FALSE)

infer_obj <- infercnv::run(
  infer_obj,
  cutoff            = 0.1,    # 文章指定
  out_dir           = out_dir,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = FALSE,  # 文中未使用 HMM
  noise_filter      = 0       # 文章指定
)
```

# 4️⃣ 手动实现文章算法：截断 + 101 滑窗 + CNV score
4.1 读取 inferCNV 输出 CNV 矩阵
```r
cnv_file <- file.path(out_dir, "infercnv.observations.txt")
cnv_mat <- read.table(cnv_file,
                      header = TRUE, row.names = 1,
                      check.names = FALSE)
cnv_mat <- as.matrix(cnv_mat)
```
4.2 表达截断到 [-3, 3]
```
cnv_trunc <- cnv_mat
cnv_trunc[cnv_trunc >  3] <-  3
cnv_trunc[cnv_trunc < -3] <- -3
```
4.3 按染色体排序并进行 101 基因滑动窗口平滑
```r
gene_order <- read.table(gene_order_file, header = TRUE, sep = "\t")
colnames(gene_order)[1:4] <- c("gene", "chr", "start", "end")

gene_order_use <- gene_order[gene_order$gene %in% rownames(cnv_trunc), ]
gene_order_use <- gene_order_use[order(gene_order_use$chr, gene_order_use$start), ]

cnv_trunc <- cnv_trunc[gene_order_use$gene, ]

window_size <- 101

chr_list <- split(gene_order_use, gene_order_use$chr)

smooth_list <- lapply(chr_list, function(df_chr) {
  genes_chr <- df_chr$gene
  mat_chr   <- cnv_trunc[genes_chr, , drop = FALSE]
  
  smoothed_chr <- apply(mat_chr, 2, function(x) {
    zoo::rollapply(x, width = window_size,
                   FUN = mean, align = "center", fill = NA)
  })

  rownames(smoothed_chr) <- genes_chr
  smoothed_chr
})

cnv_smooth <- do.call(rbind, smooth_list)
cnv_smooth <- cnv_smooth[complete.cases(cnv_smooth), ]
```
4.4 计算 CNV score（平方和）
```
CNV_score <- colSums(cnv_smooth^2)

cnv_score_df <- data.frame(
  cell      = colnames(cnv_smooth),
  CNV_score = CNV_score
)

# 写入 Seurat meta
seu$CNV_score <- cnv_score_df$CNV_score[match(colnames(seu), cnv_score_df$cell)]
```
# 5️⃣ 结果可视化（示例）
```r
library(ggplot2)

meta_epi_sub <- meta_epi[cells_use, ]
meta_epi_sub$CNV_score <- cnv_score_df$CNV_score[match(rownames(meta_epi_sub), cnv_score_df$cell)]
meta_epi_sub$group_infer <- ifelse(rownames(meta_epi_sub) %in% ref_cells,
                                   "Baseline_EN", "ET_query")

ggplot(meta_epi_sub, aes(x = group_infer, y = CNV_score, fill = group_infer)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_bw()
```
# 🔔 注意事项（务必阅读）

✔ 1. 表达矩阵格式必须是 gene × cell

inferCNV 要求：
	•	行名 = 基因
	•	列名 = 细胞条形码
	•	注释文件 cell 名必须完全一致（包括“-1” 等 10x 后缀）

⸻

✔ 2. gene_order_file 必须与你的表达基因一致

如果匹配不到基因：
	•	inferCNV 会丢掉绝大部分行
	•	CNV heatmap 会看起来“空空的”
	•	CNV score 失真

务必确保统一 hg19/hg38、gene symbol/Ensembl ID。

⸻

✔ 3. baseline（参考组）务必生物合理

文章逻辑：
	•	癌旁 EN 上皮 → 无 CNV → 可作为 baseline
	•	20% 随机抽取是为了避免 batch effect & 速度慢

如果 baseline 选不好，malignant vs non-malignant 会彻底乱掉。

⸻

✔ 4. 第一次运行建议关掉 HMM
	•	HMM 会非常耗时
	•	容易出错
	•	debug 更麻烦

跑通流程后再考虑开启。

⸻

✔ 5. 内存占用较大（单细胞数量越多越明显）

如果 Seurat 有 50k+ 上皮细胞，请考虑：
	•	先 subsample 部分用于测试
	•	再跑全量

⸻

✔ 6. inferCNV 输出文件名可能因版本不同略有改变

需要手动确认：
	•	infercnv.observations.txt
	•	或 infercnv.observations_denoised.txt

然后在脚本中对应修改。
