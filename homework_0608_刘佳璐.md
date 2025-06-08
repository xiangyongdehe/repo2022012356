### Machine Learning with R 刘佳璐
#### 8.1
```r
##### 安装必要包
library(glmnet)    # 正则化逻辑回归
library(caret)     # 数据划分、调参、特征选择等
library(pROC)      # ROC曲线绘图与AUC计算
library(ggplot2)   # PCA可视化
```

##### 数据读取和预处理
``` r
# 读取CSV
data <- read.csv("/Users/xiangyongdehe/Desktop/qPCR_data.csv")

# 标签和特征提取
y <- factor(data[[13]], levels = c("NC", "HCC"))
x <- data[, 2:12]  # 11个基因表达量

# 将所有列转为数值型
x <- apply(x, 2, as.numeric)

# 用均值填补缺失值
feature.mean <- colMeans(x, na.rm = TRUE)
x[is.na(x)] <- matrix(rep(feature.mean, each = length(y)), nrow = length(y))[is.na(x)]

# Z-score标准化
x <- scale(x, center = TRUE, scale = TRUE)
```
##### PCA可视化
```r
pca <- prcomp(x)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Label = y)

ggplot(pca_df, aes(PC1, PC2, color = Label)) +
  geom_point(size = 2.5) +
  theme_minimal() +
  labs(title = "PCA Visualization of qPCR Data")
```
PCA可视化的结果为
![image](https://github.com/user-attachments/assets/d922b0c8-5659-4bc1-a665-ab5043a4e3a2)

##### 数据集划分（80%训练，20%测试）
```r
rfFuncs$summary <- twoClassSummary
rfe_ctrl <- rfeControl(functions = rfFuncs, method = "boot", number = 10, verbose = TRUE)

set.seed(666)
rfe.results <- rfe(x.train, y.train,
                   sizes = 2:11,
                   rfeControl = rfe_ctrl,
                   metric = "ROC")

selected.features <- predictors(rfe.results)
print(selected.features)
```
##### 特征选择 RFE 基于随机森林
```r
rfFuncs$summary <- twoClassSummary
rfe_ctrl <- rfeControl(functions = rfFuncs, method = "boot", number = 10, verbose = TRUE)

set.seed(666)
rfe.results <- rfe(x.train, y.train,
                   sizes = 2:11,
                   rfeControl = rfe_ctrl,
                   metric = "ROC")

selected.features <- predictors(rfe.results)
print(selected.features)  # 输出选中的基因名
```
输出结果为
```
> print(selected.features)  # 输出选中的基因名
[1] "miR.122"          "SNORD3B"          "hsa_circ_0073052"
[4] "LINC01226"        "HULC"
```
##### 模型训练+调参数
```r
params.grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = c(0, 0.01, 0.1, 1))

tr.ctrl <- trainControl(method = "cv",
                        number = 5,
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE)

set.seed(666)
cv.fitted <- train(x.train[, selected.features],
                   y.train,
                   method = "glmnet",
                   family = "binomial",
                   metric = "ROC",
                   tuneGrid = params.grid,
                   preProcess = NULL,
                   trControl = tr.ctrl)

print(cv.fitted$bestTune)
```
##### 模型评估+ROC曲线绘制
```r
# 在测试集上预测概率
y.test.prob <- predict(cv.fitted, newdata = x.test[, selected.features], type = "prob")

# 生成 ROC 对象并绘图
roc.curve <- roc(y.test, y.test.prob[, "HCC"])
plot(roc.curve, print.auc = TRUE, main = "ROC Curve for qPCR HCC Classifier")
```
ROC曲线：
![image](https://github.com/user-attachments/assets/43c18317-129a-4763-bdbb-97f9a00f3d9c)

#### 8.2
##### 随机森林中树的数量是不是一个需要通过交叉验证调整的超参数?为什么?
在随机森林模型中，树的数量（ntree）确实是一个可以调整的超参数。它决定了整个森林中包含多少棵决策树，从而影响模型的稳定性和预测性能。当树的数量较少时，模型的预测结果可能不稳定，受样本划分影响较大；而当树的数量增加时，模型的结果会趋于稳定，方差降低，泛化能力增强。不过，树的数量并不是越多越好，因为树太多会显著增加计算时间，带来不必要的资源消耗。

在实际应用中，ntree 一般不作为最优先的调参对象。通常的做法是设定一个足够大的固定值，例如 500 或 1000，使模型达到性能收敛，再将调参重点放在影响更大的参数上，比如 mtry（每次分裂时随机选择的特征数量）以及树的最大深度等。因此，虽然 ntree 是一个可以调节的参数，但在实践中往往作为一个辅助性参数设定，而不需要通过交叉验证精细搜索。

##### 请问什么是随机森林的out-of-bag (OOB) error?它和bootstrapping有什么关系?
随机森林的 out-of-bag（OOB）误差是一种用于评估模型性能的无偏估计方法。它的原理与随机森林内部的 **bootstrapping（自助抽样）**密切相关。在构建每棵树时，随机森林会从原始训练集有放回地抽取样本来训练单棵决策树，大约有 63% 的样本被抽中用作训练，而剩下约 37% 没被抽中的样本就称为 out-of-bag 样本。对于每一棵树而言，其对应的 OOB 样本可视为“未见数据”，因此可以用来测试这棵树的预测能力。

当整个随机森林训练完成后，对于训练集中的每个样本，我们可以找到所有未在该样本参与训练的树，即该样本作为 OOB 的那部分树，并用这部分树来对它进行预测，然后将所有预测结果投票或求平均，得到该样本的 OOB 预测值。将所有样本的 OOB 预测与真实标签进行比较，最终可以计算出整个模型的 OOB误差率（OOB error），它反映了模型的泛化能力，通常与使用交叉验证得到的误差非常接近。

因此，OOB error 的本质是一种由 bootstrapping 自动带来的“内置交叉验证”，它无需额外划分验证集，也无需显式执行 K 折交叉验证，非常适合评估随机森林这类集成模型的性能。





