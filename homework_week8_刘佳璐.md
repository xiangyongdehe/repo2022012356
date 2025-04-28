## homework 8 刘佳璐
### 一、KO分析
#### 1.从wt.light.vs.dark.all.txt(这是我们在差异表达一节获得的野生型的结果)中选取显著上调的(FDR<0.05, logFC>1)的基因进行GO分析。
```r
# 读取数据
data <- read.table("/Users/xiangyongdehe/Desktop/wt.light.vs.dark.all.txt", header=TRUE, sep="\t")

# 筛选显著上调基因(FDR < 0.05且logFC > 1)
sig_up_genes <- subset(data, FDR < 0.05 & logFC > 1)

# 提取基因ID列表
gene_list <- rownames(sig_up_genes)

# 检查是否已安装
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  print("clusterProfiler 未安装")
} else {
  print("clusterProfiler 已安装")
}

# 准备输入数据
sig_genes <- rownames(subset(data, padj < 0.05 & log2FoldChange > 1))

# GO富集分析
ego <- enrichGO(gene = sig_genes,
                universe = rownames(data),
                OrgDb = org.At.tair.db,
                keyType = "TAIR",  # 直接使用TAIR ID
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
# 保存结果到文件write.csv(as.data.frame(ego), "GO_enrichment_results.csv", row.names = FALSE)
# 查看结果
head(ego)
dotplot(ego)
```
##### 结果csv文件在压缩包里

##### 输出的dotplot图![image](https://github.com/user-attachments/assets/e10a75f5-7883-414b-9fb0-da53f62c9c9f)
##### 输出的barplot图 ![image](https://github.com/user-attachments/assets/b263cabd-6567-44df-af25-31674f5f39ae)

#### 2.请问上面的例子中， Fold Enrichment和P value是如何计算的? 请写出公式，并解释原理。此外，在定义显著富集的 GO terms 时为什么一般不是参考P value的大小，而是要计算一个 FDR来做为参考？

##### Fold Enrichment 计算公式与原理
<img width="445" alt="Screenshot 2025-04-28 at 20 20 45" src="https://github.com/user-attachments/assets/971ffb05-326e-4ca7-bcbf-c259bcd6491d" />
原理：​GeneRatio (k/n): 差异基因中属于该 GO term 的比例。
​BgRatio (K/N): 基因组中所有基因属于该 GO term 的背景比例。
Fold Enrichment > 1 表示该 GO term 在差异基因中富集（比例高于背景）。

##### P value 计算公式与原理
GO 富集通常使用 ​超几何检验​（Hypergeometric Test），其 P value 计算公式为：
<img width="225" alt="Screenshot 2025-04-28 at 20 21 33" src="https://github.com/user-attachments/assets/c8a312f3-8b8b-45ae-b91b-e5f51c03541a" />
原理是检验差异基因在 GO term 中的富集是否显著高于随机抽样预期，P value < 0.05 表示富集显著（但需校正多重检验）。

##### 为什么用 FDR 而非 P value？ 
GO 分析同时测试数百至数千个 term，导致假阳性率飙升。例如，测试 1000 个 term 时，即使所有 term 都不显著，仍有约 50 个 term 的 P value < 0.05（1000×0.05）。  
而FDR为期望的假阳性结果比例，比直接使用 P value 更保守，避免过度解读随机噪声。

**总结**
| ​**指标**         | ​**计算方式**              | ​**用途**                          |  
|------------------|--------------------------|-----------------------------------|  
| Fold Enrichment  | GeneRatio / BgRatio       | 量化富集程度（倍数变化）            |  
| P value          | 超几何检验                 | 原始显著性（未校正）                |  
| FDR              | Benjamini-Hochberg 校正   | 控制多重检验假阳性，推荐阈值 < 0.05  |  

**实际分析中**：  
- 优先参考 ​**FDR**​（如 `p.adjust` 列），其次才是 P value。  
- Fold Enrichment 辅助解释生物学意义（如 Fold > 2 且 FDR < 0.05 的 term 可信度高）。

### 二、KO分析

```r
# 读取数据
data <- read.delim("/Users/xiangyongdehe/Desktop/wt.light.vs.dark.all.txt", header=TRUE, stringsAsFactors=FALSE)

# 查看数据结构
head(data)
# 筛选显著上调基因
up_genes <- subset(data, padj < 0.05 & log2FoldChange > 1)

# 查看筛选结果
dim(up_genes)
head(up_genes)

# 提取基因ID（假设基因ID是行名）
gene_symbols <- rownames(up_genes)

# 将拟南芥基因ID转换为ENTREZID
gene_ids <- bitr(gene_symbols, 
                 fromType = "TAIR",  # 拟南芥基因ID类型
                 toType = "ENTREZID", 
                 OrgDb = "org.At.tair.db")

# 进行KEGG富集分析（拟南芥的organism代码是"ath"）
kegg_result <- enrichKEGG(
  gene = gene_symbols,       # 直接使用原始TAIR ID
  organism = "ath",
  keyType = "kegg",          # 明确指定ID类型
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 查看结果
head(kegg_result)

print("KEGG结果摘要：")
head(kegg_result)

# 保存结果
write.csv(as.data.frame(kegg_result), "KEGG_enrichment_results.csv", row.names=FALSE)

# 结果可视化
# 绘制条形图
barplot(kegg_result, showCategory=20)
# 绘制点图
dotplot(kegg_result, showCategory=20)
# 绘制基因-通路网络图
cnetplot(kegg_result, categorySize="pvalue", foldChange=up_genes$logFC)

```
