
# 生物信息学学习记录 | 刘佳璐

# 第一课课堂笔记
# Question/Hypothesis-driven Science
## 4 steps of Bioinformatics
### step 0 Question: Biological/Medical Knowledge
### step 1 Information: Biological/Medical Data
### step 2 Analysis: Data clean & Feature Extraction
* **NGS Data Analysis**  
  Variation | Sequencing Methods | Bioinformatics Tools
  --- | --- | ---
  Abundance | RNA-seq | DEGseq2,EdgeR,Cufflinks
  Splicing | RNA-seq | rMATs,TOPHAT/Cufflinks
  APA | RNA-seq/PAT-seq | DaPars,APAtrap
  Translation | Ribo-seq | RiboWave,RiboTaper,ORFscore
  Degradation | Degradome-Seq,cfRNA-seq | sPARTA,dPeak
  Modification | m6A-seq,MeRIP-seq,miCLIP | m6aViewer,MeRIP-PF
  Structure | icSHAPE,SHAPE-map,DMS-seq | RNAstructure,RNAfold,RME
  RNA-protein interaction | HITS-CLIP,PAR-CLIP,iCLIP,eCLIP | Piranha,PARalyzer,CIMS
### step 3 Modeling: Probabilistic Model & Computational Algorithm
* **Probabilistic Model:**  
  A mathematical framework that uses **probability theory** to represent uncertainty in systems or processes. It describes relationships between random variables through probability distributions (e.g., Bayesian networks, hidden Markov models).  
* **computational Algorithm:**  
  A **step-by-step procedure** to solve a computational problem or perform a task. It defines a finite sequence of operations to transform inputs into outputs (e.g., sorting algorithms, dynamic programming).
* **Comparation**
  Aspect | Probabilistic Model | Computational Algorithm
  --- | --- | ---
  Primary Goal | Model uncertainty and stochastic relationships | Solve problems efficiently
  Core Components | Random variables, probability distributions | Steps, loops, conditionals
  Evaluation | Metric Likelihood, posterior predictive accuracy |Time/space complexity, correctness
  Output | Type Probability distributions or confidence scores | Deterministic or approximate solutions
  Key Applications | Risk assessment, gene prediction, NLP | Data sorting, optimization, simulations
  Mathematical Basis | Probability theory, statistics | Discrete mathematics, graph theory
  Implementation Focus | Parameter estimation, inference methods | Code optimization, resource management

  
# Big Data-driven Science
## 4 steps of Bioinformatics
### step 0 Information: Biological/Medical Data
### step 1 Analysis: Data clean & Feature Extraction
### step 2 Modeling: Probabilistic Model & Computational Algorithm
### step 3 Question: Biological/Medical Knowledge

# 本学期生信学习计划
## 代码基础部分
* week 1-4: Linux
* week 5-16: R
* week 11-16: python
* 每天固定学习至少两个小时
## 实操部分
* NGS data analysis
* 自己从seurat网站上学习组学分析，各种测序分析画图的基本流程
* 每两周读一篇生信文献，根据文献数据复现分析结果图
