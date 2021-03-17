# 气体特性模拟

## 模拟物理参数
- 漂移速度 (Drift Velocity)
- 横向扩散、纵向扩散 (Diffusion Coefficient)
- 汤生系数 (Townsend Coefficient)
- 吸附系数 (Attachment Coefficient)

## 程序说明
程序包含三个文件：
- generate.C
- read.C
- makefile
### 生成气体文件
气体文件的后缀为 `.gas`，创建 `result` 与 `runlog` 文件夹，`result` 用来存放 `gas` 文件，`runlog` 用来存放运行日志。运行日志包含所有的计算结果。
```bash
make generate
./generate [-g1 Gas1] [-f1 Fraction1] [-g2 Gas2] [-f2 Fraction2] [-g3 Gas3] [-f3 Fraction3] [-p Pressure] [-t Temperature]
```
最多可以设置三种气体，可以设置不同比分、不同气压 (atm)、不同温度 (degree centigrade)。

### 结果绘图

#### 读取 `gas` 文件

```bash
make read
./read par
```
`par` 为 `gas` 文件的路径。

#### 读取 `log` 文件

通过 `shell` 脚本从 `log` 文件中提取信息到 `txt` 文件，再通过 `ROOT` 程序绘图。
```bash
#!/bin/bash
filename=xxx.log
EF=`grep -nr "ELECTRIC FIELD" ./runlog/$filename | awk -F ' ' '{print $5}'`
DV=`grep -nr "along E:" ./runlog/$filename | awk -F ' ' '{print $6}'`
TD=`grep -nr "Transverse diffusion" ./runlog/$filename | awk -F ' ' '{print $4}'`
LD=`grep -nr "Longitudinal diffusion" ./runlog/$filename | awk -F ' ' '{print $4}'`
TC=`grep -nr "Townsend coefficient:" ./runlog/$filename | awk -F ' ' '{print $4}'`
AC=`grep -nr "Attachment coefficient:" ./runlog/$filename | awk -F ' ' '{print $4}'`
echo $EF $DV $TD $LD $TC $AC > ./result/data.txt
```

