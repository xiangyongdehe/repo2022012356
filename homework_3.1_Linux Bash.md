# 1.3 Linux Bash 学习笔记
##  example 1: run a bash script
```bash
touch bash.sh #新建文件
chmod u+x bash.sh #    为该文件添加可执行权限
#编辑bash.sh，内容如下：
vim bash.shv:w
#按下 i 键进入插入模式（Insert mode），此时你可以开始编辑文件。
#输入以下内容：
#!/bin/bash
echo "hello bash"
exit 0
#按下 Esc 键退出插入模式。
#输入 :wq 并按下 Enter 键，保存文件并退出 vim。
# 执行脚本。这种方法需要脚本文件有可执行的权限，这是我们前面执行命令chmod u+x bash.sh的目的所在
./bash.sh
```
## if 语句块
```bash
#!/bin/bash

# 提示用户输入一个值
echo -n "please input a number:"

# 保存用户输入的值到num中
read num

if [ "$num" -lt "0" ];then
# 小于0,则输出“negtive number” -lt 小于 -eq 等于 -gt大于
echo "negtive number"
elif [ "$num" -gt "0" ];then
# 大于0,则输出“positive number”
echo "positive number"
else
# 大于0,则输出"number zero"
echo "number zero"
fi
exit 0
```


