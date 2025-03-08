#  homework 
## 参考和学习本章内容，写出一个 bash 脚本，可以使它自动读取一个文件夹（例如 bash_homework/）的内容，将该文件夹下文件的名字输出到 filenames.txt, 子文件夹的名字输出到 dirname.txt 。

## 编写脚本
```bash
xiangyongdehe@xiangyoehedeAir ~ % docker run -it -v ~/Desktop/bash_homework:/home/test/linux/bash_homework xfliu1995/bioinfo_tsinghua:2 /bin/bash
WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested
test@aca2dba3acc2:~$ cd /home/test/linux
test@aca2dba3acc2:~/linux$ ls
1.gtf.gz  bash_homework  file
#挂载bash_homework到容器的/home/test/linux

touch file.sh
chmod u+x file.sh
vim file.sh
#随后输入i编写脚本

# 定义目标文件夹
TARGET_DIR="bash_homework"

# 检查目标文件夹是否存在
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory $TARGET_DIR does not exist."
    exit 1
fi

# 清空或创建输出文件
> filenames.txt
> dirnames.txt

# 遍历目标文件夹
for item in "$TARGET_DIR"/*; do
    if [ -f "$item" ]; then
        # 如果是文件，将文件名写入 filenames.txt
        echo "$(basename "$item")" >> filenames.txt
    elif [ -d "$item" ]; then
        # 如果是子文件夹，将文件夹名写入 dirnames.txt
        echo "$(basename "$item")" >> dirnames.txt
    fi
done

echo "File names have been saved to filenames.txt."
echo "Directory names have been saved to dirnames.txt."

#按esc，并且输入：wq保存脚本
```
## 运行脚本
```bash
./file.sh
File names have been saved to filenames.txt.
Directory names have been saved to dirnames.txt.
```
## 检查输出结果
```bash
#查看文件
cat filenames.txt
#输出
a1.txt
a.txt
b1.txt
bam_wig.sh
b.filter_random.pl
c1.txt
chrom.size
c.txt
d1.txt
dir.txt
e1.txt
f1.txt
human_geneExp.txt
if.sh
image
insitiue.txt
mouse_geneExp.txt
name.txt
number.sh
out.bw
random.sh
read.sh
test3.sh
test4.sh
test.sh
test.txt
wigToBigWig

#查看文件
cat dirnames.txt
# 输出
a-docker
app
backup
bin
biosoft
c1-RBPanno
datatable
db
download
e-annotation
exRNA
genome
git
highcharts
home
hub29
ibme
l-lwl
map2
mljs
module
mogproject
node_modules
perl5
postar2
postar_app
postar.docker
RBP_map
rout
script
script_backup
software
tcga
test
tmp
tmp_script
var
x-rbp
