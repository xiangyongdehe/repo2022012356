docker run -v /Users/xiangyongdehe/Desktop/test_command.gtf:/data/test_command.gtf xfliu1995/bioinfo_tsinghua:2 # 运行容器
docker run -it -v /Users/xiangyongdehe/Desktop/test_command.gtf:/data/test_command.gtf xfliu1995/bioinfo_tsinghua:2 /bin/bash # 验证文件是否挂载
test@fa074a67909f:~$ cat /data/test_command.gtf # 验证文件是否能打开查看
chr_IV	ensembl	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1";
chr_IV	ensembl	transcript	802	2953	.	+	.	gene_id "YDL248W"; gene_version "1";
chromosome_IV	ensembl	exon	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1";
chromosome_IV	ensembl	CDS	1802	950	.	+	0	gene_id "YDL248W"; gene_version "1";
chr_IV	ensembl	start_codon	1802	1804	.	+	0	gene_id "YDL248W"; gene_version "1";
chromosome_IV	ensembl	stop_codon	2951	2953	.	+	0	gene_id "YDL248W"; gene_version "1";
chromosome_IV	ensembl	gene	762	3836	.	+	.	gene_id "YDL247W-A"; gene_version "1";
chr_IV	ensembl	transcript	3762	836	.	+	.	gene_id "YDL247W-A"; gene_version "1"; # 输出行

##task 1
```bash
$ wc -l /data/test_command.gtf  # 输入：统计文件行数
8 /data/test_command.gtf        # 输出：行数结果
$ wc -m /data/test_command.gtf  # 输入：统计文件字符数
636 /data/test_command.gtf      # 输出：字符数结果

##task 2
```bash
~$ grep -E '^chr_.*\bYDL248W\b' test_command.gtf # 输入：筛选并输出示例文件中以 chr_ 起始，并且基因id为 YDL248W 的行
grep: test_command.gtf: No such file or directory
test@fa074a67909f:~$ grep -E '^chr_.*YDL248W' /data/test_command.gtf
chr_IV	ensembl	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1";
chr_IV	ensembl	transcript	802	2953	.	+	.	gene_id "YDL248W"; gene_version "1";
chr_IV	ensembl	start_codon	1802	1804	.	+	0	gene_id "YDL248W"; gene_version "1"; # 输出：筛选结果

##task 3
```bash
~$ sed 's/chr_/chromosome_/g' /data/test_command.gtf | cut -f 1,3,4,5 # 将示例文件中的 chr_ 替换为 chromosome_ 并输出每行的第1，3，4，5列
chromosome_IV	gene	1802	2953
chromosome_IV	transcript	802	2953
chromosome_IV	exon	1802	2953
chromosome_IV	CDS	1802	950
chromosome_IV	start_codon	1802	1804
chromosome_IV	stop_codon	2951	2953
chromosome_IV	gene	762	3836
chromosome_IV	transcript	3762	836 # 输出替换结果

##task 4
```bash
awk 'BEGIN {OFS="\t"} {tmp=$2; $2=$3; $3=tmp; print}' /data/test_command.gtf | sort -k4,4n -k5,5n # 通过man命令以及更多的资料学习简单的 awk 命令，尝试互换示例文件的第2列和第3列，并且对输出结果利用 sort 命令依照第4和第5列数字大小排序，将最终结果输出到result.gtf文件中。
chromosome_IV	gene	ensembl	762	3836	.	+	.	gene_id"YDL247W-A";	gene_version	"1";
chr_IV	transcript	ensembl	802	2953	.	+	.	gene_id"YDL248W";	gene_version	"1";
chromosome_IV	CDS	ensembl	1802	950	.	+	0	gene_id"YDL248W";	gene_version	"1";
chr_IV	start_codon	ensembl	1802	1804	.	+	0	gene_id"YDL248W";	gene_version	"1";
chr_IV	gene	ensembl	1802	2953	.	+	.	gene_id	"YDL248W";	gene_version	"1";
chromosome_IV	exon	ensembl	1802	2953	.	+	.	gene_id"YDL248W";	gene_version	"1";
chromosome_IV	stop_codon	ensembl	2951	2953	.	+	0	gene_id	"YDL248W";	gene_version	"1";
chr_IV	transcript	ensembl	3762	836	.	+	.	gene_id"YDL247W-A";	gene_version	"1"; # 输出结果
awk 'BEGIN {OFS="\t"} {tmp=$2; $2=$3; $3=tmp; print}' /data/test_command.gtf | sort -k4,4n -k5,5n > result.gtf # 将最终结果输出到result.gtf文件中

##task 5
~$ ls -l /data/test_command.gtf # 查看文件的当前权限
-rw-r--r-- 1 test test 636 Mar  4 07:07 /data/test_command.gtf # 输出：文件所有者：读、写 (rw-)；用户组：读 (r--)；其他用户：读 (r--)
chmod 754 /data/test_command.gtf
ls -l /data/test_command.gtf # 修改文件权限
-rwxr-xr-- 1 test test 636 Mar  4 07:07 /data/test_command.gtf# 文件所有者（test）现在可以读、写、执行文件，用户组（test）可以读和执行文件，而其他用户只能读取文件
