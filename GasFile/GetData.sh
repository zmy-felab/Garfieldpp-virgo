#!/bin/bash
# 从log文件中提取数据，只适用于3.0及master版本输入的log文件。
# sh GetData.sh xxx，xxx为log文件名，不带后缀。

# 提取数据信息
EF=`grep -nr "ELECTRIC FIELD" ./runlog/$1.log | awk -F '=' '{print $2}' | awk '{print $1}'`
DV=`grep -nr "along E:" ./runlog/$1.log | awk -F ' ' '{print $6}'`
TD=`grep -nr "Transverse diffusion" ./runlog/$1.log | awk -F ' ' '{print $4}'`
LD=`grep -nr "Longitudinal diffusion" ./runlog/$1.log | awk -F ' ' '{print $4}'`
TC=`grep -nr "Townsend coefficient:" ./runlog/$1.log | awk -F ' ' '{print $4}'`
AC=`grep -nr "Attachment coefficient:" ./runlog/$1.log | awk -F ' ' '{print $4}'`

# 判断runtxt文件夹是否存在
if [ ! -d "./runtxt" ]
then
    mkdir runtxt
fi

# 保存数据，每个结果一行
echo $EF > ./runtxt/temp
echo $DV >> ./runtxt/temp
echo $TD >> ./runtxt/temp
echo $LD >> ./runtxt/temp
echo $TC >> ./runtxt/temp
echo $AC >> ./runtxt/temp

# 行列互换
cat ./runtxt/temp | awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(i=1;i<=NF;i++){for(j=1;j<=NR;j++)printf a[j,i]" ";print ""}}' > ./runtxt/$1.txt

# 删除临时文件
rm -rf ./runtxt/temp
