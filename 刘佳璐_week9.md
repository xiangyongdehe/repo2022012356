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
docker exec -it bioinfo_tsinghua bash
cd /home/test/chip-seq/

##### step 2 创建tag目录
makeTagDirectory ip_tagdir ip.chrom_part.bam
makeTagDirectory input_tagdir input.chrom_part.bam

##### step 3 peak calling
findPeaks ip_tagdir -style factor -i input_tagdir -o peaks0.txt

##### step 4 motif分析
findMotifsGenome.pl peaks0.txt sacCer2 motif_output0 -len 8

findMotifsGenome.pl peaks0.txt sacCer2 motif_output_adjusted \
  -len 10 \        
  -size 200 \       
  -p 0.1 \          
  -mask    
  ==============
 findMotifsGenome.pl peaks0.txt sacCer2 motif_output0 -len 10

  ls -lh motif_output0 \

# 从容器复制到当前目录（在容器内执行）
docker cp c18811589b5b:~/chip-seq/homework/motif_output0/ ./motif_output_local/

# 输出文件是

