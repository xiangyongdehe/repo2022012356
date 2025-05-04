#### 请解释在ChIP-seq实验中为什么一般都要平行做一个 control （通常叫 input）的实验。
#### 请解释 findPeaks 和 findMotifsGenome.pl 主要参数的含义。
#### 我们在容器的/home/test/chip-seq/homework目录中提供了酵母Snf1蛋白CHIP-seq的bam文件，ip.chrom_part.bam为IP实验数据，input.chrom\_part.bam为背景数据。请大家从这两个文件出发，用homer重复本章中介绍的peak calling和motif finding分析。请大家提交找到的motif的截图，以及Fold Change (vs Control) >=8且p-value (vs Control) < 10的−8次方的peaks(建议放在同一个文件中提交)。
##### step 1 安装docker 并挂载容器
```bash
docker load -i "/Users/xiangyongdehe/Desktop/bioinfo/bioinfo_motif.2.0.tar"
open /Users/xiangyongdehe/Desktop/bioinfo folder/bioinfo_clip_seq.tar: no such file or directory
e00b027e565f: Loading layer  2.599GB/2.599GB
# 运行容器并挂载桌面上bioinfo文件夹里的文件到容器内
docker run -it \
  -v /Users/xiangyongdehe/Desktop/bioinfo/:/home/test/chip-seq/homework \
  zwhbio2017/clip-seq /bin/bash
```
##### step 2 peaking calling 分析
```bash
# 创建HOMER tag目录 便于分析
makeTagDirectory ip_tagdir/ ip.chrom_part.bam
makeTagDirectory input_tagdir/ input.chrom_part.bam

# 安装 sacCer3 基因组
perl /usr/local/bin/HOMER/configureHomer.pl -install sacCer3

# 验证安装
ls /usr/local/bin/HOMER/data/genomes/sacCer3/

findPeaks ip.chrom_part.bam \
  -style factor \
  -i input.chrom_part.bam \
  -o peaks.txt \
  -log2FoldThresh 3 \  # 对应Fold Change ≥ 8
  -pvalueThresh 1e-8   # p-value阈值
```

##### step 3 提取显著peak
```bash
# 提取FC≥8且p-value<1e-8的peaks
awk 'BEGIN{OFS="\t"; print "PeakID\tChr\tStart\tEnd\tStrand\tFoldChange\tp-value"}
     NR>1 && $6>=8 && $8<1e-8 {print $1,$2,$3,$4,$5,$6,$8}' peaks.txt > significant_peaks.txt

# 统计显著peaks数量
wc -l significant_peaks.txt
```

##### step 4 motif 分析
```bash
findMotifsGenome.pl significant_peaks.txt sacCer3 motif_output \
  -size 100 \          # 分析peaks中心±100bp
  -mask \              # 屏蔽重复序列
  -bg input.chrom_part.bam  # 使用input作为背景
```

##### step 5 结果处理与可视化
```bash
# 1. 获取top motif截图
convert motif_output/homerResults/logo1.png -resize 50% top_motif.png

# 2. 生成结果报告
echo "## 酵母Snf1蛋白ChIP-seq分析结果" > results.md
echo "### 显著peaks统计" >> results.md
awk 'END{print "共发现"NR-1"个显著peaks(FC≥8,p<1e-8)"}' significant_peaks.txt >> results.md
echo "### Top Motif" >> results.md
echo "![Top Motif](top_motif.png)" >> results.md

# 3. 打包结果文件
zip submission.zip significant_peaks.txt top_motif.png results.md
```
findPeaks ip.chrom_part.bam \
>   -style factor \
>   -i /home/test/chip-seq/homework/input_tagdir/ \
>   -o /home/test/chip-seq/homework/peaks.txt \
>   -log2FoldThresh 3 \
>   -pvalueThresh 1e-8
bash: /usr/local/bin/HOMER/bin/findPeaks: No such file or directory

