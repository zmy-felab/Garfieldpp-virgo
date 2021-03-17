# 气体特性模拟

## 模拟物理参数
- 漂移速度 (Drift Velocity)
- 横向扩散、纵向扩散 (Diffusion Coefficient)
- 汤生系数 (Townsend Coefficient)
- 吸附系数 (Attachment Coefficient)

## 程序说明

### 生成气体文件
气体文件的后缀为 `.gas`，创建 `result` 与 `runlog` 文件夹，`result` 用来存放 `gas` 文件，`runlog` 用来存放运行日志。运行日志包含所有的运行结果。
```bash
make generate
./generate [-g1 Gas1] [-f1 Fraction1] [-g2 Gas2] [-f2 Fraction2] [-g3 Gas3] [-f3 Fraction3] [-p Pressure] [-t Temperature]
```
最多可以设置三种气体，不同比分，不同气压 (atm)，不同温度 (degree centigrade)。

### 物理量

- 读取气体文件
    ```bash
    make read
    ./read par
    ```
    `par` 为 `gas` 文件的路径。

- 读取 `log` 文件
    通过 `shell` 命令从 `log` 文件中提取信息，再通过 `ROOT` 绘图。
    ```bash
    grep -nr "ELECTRIC FIELD" ./runlog/xxx.log | awk -F ' ' '{print $5}' > ./result/ElectricField.txt
    grep -nr "along E:" ./runlog/xxx.log | awk -F ' ' '{print $6}' > ./result/DriftVelocity.txt
    grep -nr "Transverse diffusion" ./runlog/xxx.log | awk -F ' ' '{print $4}' > ./result/TransverseDiffusion.txt
    grep -nr "Longitudinal diffusion" ./runlog/xxx.log | awk -F ' ' '{print $4}' > ./result/LongitudinalDiffusion.txt
    grep -nr "Townsend coefficient:" ./runlog/xxx.log | awk -F ' ' '{print $4}' > ./result/TownsendCoefficient.txt
    grep -nr "Attachment coefficient:" ./runlog/xxx.log | awk -F ' ' '{print $4}' > ./result/AttachmentCoefficient.txt
    ```

