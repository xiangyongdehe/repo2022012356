```r
library(Seurat)
library(dplyr)
library(patchwork)

# 1) 读取10x数据（cellranger 输出目录）
pbmc.data <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix/")  # 或 Read10X_h5("matrix.h5")

# 2) 创建对象
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# 3) 计算线粒体比例（人类基因通常用 ^MT- ，小鼠用 ^mt-）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4) 质控与筛选（阈值可按项目微调）
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 5A) 经典规范化工作流（LogNormalize + HVGs + ScaleData）
pbmc <- NormalizeData(pbmc)                                   # LogNormalize, scale.factor=1e4
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)                                        # 默认只会缩放可变基因
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
ElbowPlot(pbmc)                                                # 观察PC拐点

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# 5B) 或者直接用 SCTransform（推荐，更稳健；会替代 Normalize/FindVariable/Scale）
# pbmc <- SCTransform(pbmc, vst.flavor = "v2", verbose = FALSE)
# pbmc <- RunPCA(pbmc)
# pbmc <- RunUMAP(pbmc, dims = 1:30)
# pbmc <- FindNeighbors(pbmc, dims = 1:30)
# pbmc <- FindClusters(pbmc, resolution = 0.5)

# 6) 差异基因（建议装 presto 以提速）
# install.packages("presto")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)             # 与其他所有细胞对比
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# 7) 可视化
VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeaturePlot(pbmc, features = c("MS4A1","CD3E","GNLY","CD14","FCGR3A","LYZ","PPBP","CD8A"))

# 8) 保存
saveRDS(pbmc, file = "pbmc_tutorial_final.rds")
```
