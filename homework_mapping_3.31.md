# mapping
## 1.请阐述bowtie中利用了 BWT 的什么性质提高了运算速度？并通过哪些策略优化了对内存的需求？
Bowtie 利用 Burrows-Wheeler Transform (BWT) 的关键性质来提高运算速度，并通过多种策略优化内存需求。BWT 的核心优势在于其压缩后的高效搜索能力。它将原始序列转换为一种更易于压缩的形式，同时保留所有原始信息，使得在压缩数据上可以直接进行高效的子串搜索。Bowtie 结合 BWT 和 FM-index，利用反向搜索（Backward Search）逐步扩展查询序列的后缀匹配，避免了全序列扫描。FM-index 通过预先计算的辅助数据结构（如 rank 和 count 表）实现常数时间内的字符定位，极大减少了搜索时间。此外，BWT 的结构支持精确匹配和近似比对，动态规划结合 BWT 的跳跃机制（如双重索引策略）在保持高效的同时允许错配或插入缺失。

在内存优化方面，Bowtie 采用了几种关键策略。首先，BWT 转换后的序列通常包含长串重复字符，可通过游程编码（RLE）等压缩方法减少存储空间。其次，Bowtie 仅存储稀疏采样的后缀数组（SA），在查询时通过 LF-mapping 和局部解码恢复完整信息，大幅降低内存占用。此外，Bowtie 将大型基因组分割为多个小块，分别构建 BWT 索引，避免一次性加载整个基因组。比对时仅加载当前查询所需的索引部分，进一步减少内存峰值使用。Bowtie 还使用紧凑的数据结构，如压缩的 rank 和 count 表，利用位操作（如 bit-packed integers）存储辅助数据，而非传统整数类型。在 Bowtie 1 中，双重 BWT 索引策略（同时构建原始序列和反向互补序列的索引）通过共享部分数据结构（如 rank 表）减少冗余存储。

这些设计使得 Bowtie 在短序列比对中表现出色，既提高了运算速度，又优化了内存需求，使其成为高通量测序分析的高效工具。后续版本（如 Bowtie 2）进一步优化算法，支持更长的读长和更复杂的比对需求。

## 2.用bowtie将 THA2.fa mapping 到 BowtieIndex/YeastGenome 上，得到 THA2.sam，统计mapping到不同染色体上的reads数量(即统计每条染色体都map上了多少条reads)。
```bash
docker exec -it bioinfo_tsinghua /bin/bash
# 进入本节课的工作目录
cd /home/test/mapping

~$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
#从NCBI上下载酵母基因组
~$ gunzip GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
~$ bowtie-build GCF_000146045.2_R64_genomic.fna BowtieIndex/YeastGenome
~$ bowtie-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa BowtieIndex/YeastGenome
~$ bowtie -S BowtieIndex/YeastGenome THA2.fa THA2.sam
#使用索引进行比对
~$ echo "Chromosome mapping statistics:"
grep -v "^@"  THA2.sam | cut -f3 | sort | uniq -c | sort -nr
#输出
    189 chrIV
    168 chrXII
    124 chrVII
    101 chrXV
     79 *
     78 chrXVI
     72 chrXIII
     70 chrX
     69 chrVIII
     61 chrXIV
     54 chrXI
     51 chrII
     33 chrV
     28 chrIX
     20 chrIII
     20 chrI
     17 chrVI
     16 chrmt

```

## 3.查阅资料，回答以下问题:
### （3.1）什么是sam/bam文件中的"CIGAR string"? 它包含了什么信息?
3.1CIGAR string是 SAM/BAM 文件中的一个字段，用于描述测序读段（read）与参考基因组的比对情况。它由 <长度><操作符> 的组合构成，记录匹配（M）、插入（I）、缺失（D）、软裁剪（S）等比对细节。

### （3.2）"soft clip"的含义是什么，在CIGAR string中如何表示？
3.2Soft clip（软裁剪）​ 表示读段的部分碱基未被比对（如测序低质量区域或接头污染），但这些碱基仍保留在测序数据中。在 CIGAR string 中以 S 表示（如 2S3M 表示前 2 个碱基被软裁剪，后 3 个匹配）
### （3.3）什么是reads的mapping quality? 它反映了什么样的信息?
Mapping quality（MAPQ）​ 是 SAM/BAM 文件中一个 0–255 的数值（如 60），表示读段比对到参考基因组的可靠性。数值越高，比对越可信（如 60 通常代表唯一比对），低值（如 0）可能表示多重比对或低置信度。
### （3.4）仅根据sam/bam文件的信息，能否推断出read mapping到的区域对应的参考基因组序列? 
3.4可以推断，但需要结合 ​CIGAR string 和 ​MD tag（可选字段）​ 的信息——CIGAR string 提供比对结构（匹配/插入/缺失等），确定 read 如何比对到参考基因组。MD 标签，它记录了 read 与参考序列之间的匹配信息，包括匹配的碱基和不匹配的碱基位置。通过解析 MD 标签，可以 reconstruct 出相应的参考序列。

## 4. 请自行安装教程中未涉及的bwa软件，从UCSC Genome Browser下载Yeast (S. cerevisiae, sacCer3)基因组序列。使用bwa对Yeast基因组sacCer3.fa建立索引，并利用bwa将THA2.fa，mapping到Yeast参考基因组上，并进一步转化输出得到THA2-bwa.sam文件。
```bash
#安装bwa
git clone https://github.com/lh3/bwa.git
cd bwa; make
#从UCSC Genome Browser下载Yeast (S. cerevisiae, sacCer3)基因组序列
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
# 解压
gunzip sacCer3.fa.gz
# 重命名为 ref.fa
mv sacCer3.fa ref.fa
#使用 BWA 构建基因组索引
cd ~/bwa  
 ./bwa index ref.fa
# 运行比对 利用bwa将THA2.fa，mapping到Yeast参考基因组上，并进一步转化输出得到THA2-bwa.sam文件
./bwa mem ref.fa /Users/xiangyongdehe/./THA2.fa > THA2-bwa.sam
#查看结果文件前15行
head -n 15 THA2-bwa.sam
#输出
@SQ	SN:chrI	LN:230218
@SQ	SN:chrII	LN:813184
@SQ	SN:chrIII	LN:316620
@SQ	SN:chrIV	LN:1531933
@SQ	SN:chrIX	LN:439888
@SQ	SN:chrV	LN:576874
@SQ	SN:chrVI	LN:270161
@SQ	SN:chrVII	LN:1090940
@SQ	SN:chrVIII	LN:562643
@SQ	SN:chrX	LN:745751
@SQ	SN:chrXI	LN:666816
@SQ	SN:chrXII	LN:1078177
@SQ	SN:chrXIII	LN:924431
@SQ	SN:chrXIV	LN:784333

```
