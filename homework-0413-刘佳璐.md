# 第六次作业 刘佳璐 2022012356
```
（1）我们提供的bam文件COAD.ACTB.bam是单端测序分析的结果还是双端测序分析的结果？为什么？(提示：可以使用samtools flagstat）

（2）查阅资料回答什么叫做"secondary alignment"？并统计提供的bam文件中，有多少条记录属于"secondary alignment?" （提示：可以使用samtools view -f 获得对应secondary alignment的records进行统计）

（3）请根据hg38.ACTB.gff计算出在ACTB基因的每一条转录本中都被注释成intron的区域，以bed格式输出。并提取COAD.ACTB.bam中比对到ACTB基因intron区域的bam信息，后将bam转换为fastq文件。

提示：

写脚本把ACTB在gff中第三列为"gene"的interval放在一个bed文件中，第三列为"exon"的intervals放在另外一个bed文件中，再使用bedtools subtract。

请注意bed文件使用的是0-based coordinate，gff文件使用的是1-based coordinate。

鼓励其他实现方法，描述清楚过程即可

(4) 利用COAD.ACTB.bam计算出reads在ACTB基因对应的genomic interval上的coverage，以bedgraph格式输出。 （提示：对于真核生物转录组测序向基因组mapping得到的bam文件，bedtools genomecov有必要加-split参数。）
```

## 1）我们提供的bam文件COAD.ACTB.bam是单端测序分析的结果还是双端测序分析的结果？为什么？

```bash
test@bioinfo_docker:~/share/homework$ samtools flagstat COAD.ACTB.bam
185650 + 0 in total (QC-passed reads + QC-failed reads)
4923 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
185650 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
### paired in sequencing、read1 和 read2 是 0，这意味着bam文件是单端测序的结果。


## 2) 查阅资料回答什么叫做"secondary alignment"？并统计提供的bam文件中，有多少条记录属于"secondary alignment?" 
```bash
view -f 256 COAD.ACTB.bam | wc -l
4923
```
在二代测序（NGS）数据比对中，"Secondary Alignment"（次要对齐）是指一条测序reads（序列）比对到参考基因组上的多个位置，但只有其中一个位置被标记为主要比对（Primary Alignment）​，其余位置则被标记为次要对齐。
有4923条记录属于"secondary alignment"

## 3) 请根据hg38.ACTB.gff计算出在ACTB基因的每一条转录本中都被注释成intron的区域，以bed格式输出。并提取COAD.ACTB.bam中比对到ACTB基因intron区域的bam信息，后将bam转换为fastq文件。
```bash
awk '$3=="gene" && $9~/ACTB/{print $1"\t"$4-1"\t"$5}' hg38.ACTB.gff > ACTB_gene.bed
awk '$3=="exon"{print $1"\t"$4-1"\t"$5}' hg38.ACTB.gff > ACTB_exon.bed
bedtools subtract -a ACTB_gene.bed -b ACTB_exon.bed > ACTB_intron.bed
# 提取intron
chr7	5528185	5528280
chr7	5529982	5530523
chr7	5530627	5540675
chr7	5540771	5561851
chr7	5561949	5562389
chr7	5562828	5563713
# 输出ACTB_intron.bed的内容
bedtools intersect -a COAD.ACTB.bam -b ACTB_intron.bed -wa > COAD_ACTB_intron.bam
# 提取比对到 ACTB intron 区域的 BAM 文件
samtools fastq COAD_ACTB_intron.bam > COAD_ACTB_intron.fastq
# 将 BAM 转换为 FASTQ
```

## 4) 利用COAD.ACTB.bam计算出reads在ACTB基因对应的genomic interval上的coverage，以bedgraph格式输出。 （提示：对于真核生物转录组测序向基因组mapping得到的bam文件，bedtools genomecov有必要加-split参数。）

```bash
samtools sort COAD.ACTB.bam -o COAD_ACTB.sorted.bam
samtools index COAD_ACTB.sorted.bam
samtools view -b -L ACTB_gene.bed COAD_ACTB.sorted.bam > COAD_ACTB_region.bam
bedtools genomecov -ibam COAD_ACTB_region.bam -bga -split > ACTB_coverage.bedgraph
#利用COAD.ACTB.bam计算出reads在ACTB基因对应的genomic interval上的coverage
test@bioinfo_docker:~/share/homework$ head -5 ACTB_coverage.bedgraph 
chr7    0       5045717 0
chr7    5045717 5045731 1
chr7    5045731 5058689 0
chr7    5058689 5058695 1
chr7    5058695 5072542 0
#bedgraph的输出结果
```

## 附加作业

<img width="839" alt="Screenshot 2025-04-14 at 21 53 11" src="https://github.com/user-attachments/assets/d13949a0-7cc1-4fa1-8dd5-6ce6b77b018c" />

# 1.人类基因组的大小及组成
人类基因组大小约为32亿个碱基对（3.2Gb），其组成为：

| 区域         | 比例     | 功能说明                             |
|--------------|----------|--------------------------------------|
| 编码序列       | 约1.5%   | 负责蛋白质合成             |
| 非编码区     | 约98.5%  | 包含内含子、调控序列、非编码RNA、重复序列等 |
| 重复序列     | 约50%    | 如LINEs、SINEs、LTR等重复元素        |
| 线粒体DNA    | 约16,569 bp | 编码13种蛋白和多个tRNA/rRNA        |

![image](https://github.com/user-attachments/assets/80ec9b2d-9ae9-4e95-a026-8cc92dd2896c)

# 2. 基因中非编码RNA的最新注释与功能分类
|  RNA类型   | 基因/转录本数量（约） | 功能简述                                                                 |
|-----------|----------------------|--------------------------------------------------------------------------|
| ​**lncRNA**    | 18,000转录本         | 通过染色质重塑、miRNA海绵等机制调控基因表达，参与癌症和发育调控 |
| ​**miRNA**     | 2,000个              | 结合RISC复合体调控mRNA稳定性，介导翻译抑制或降解              |
| ​**snoRNA**    | 2,200个              | 指导rRNA/tRNA的化学修饰（如甲基化），部分参与癌症转移调控       |
| ​**snRNA**     | 1,400个              | 形成剪接体识别pre-mRNA剪接位点，U1/U2调控剪接起始              |
| ​**rRNA**      | 600个                | 构成核糖体并催化肽键形成，决定翻译保真度                     |
| ​**tRNA**      | 600个                | 转运氨基酸至核糖体，新发现通过PTMD机制介导mRNA降解          |
| ​**circRNA**   | 100,000+             | 环状结构稳定，作为miRNA海绵或蛋白互作支架                      |
| ​**piRNA**     | 30,000+              | 与PIWI蛋白结合沉默生殖细胞中转座子，维持基因组稳定性           |
| ​**eRNA**      | 50,000+              | 增强子活性调控，促进染色质环化或转录因子招募                     |

**参考数据** 
[NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/genome/gdv/)
[ENSEMBL Homo sapiens Genome](https://www.ensembl.org/Homo_sapiens/Info/Index)

