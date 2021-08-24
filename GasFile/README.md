# 气体特性模拟

## 模拟物理参数
- 漂移速度 (Drift Velocity)
- 横向扩散、纵向扩散 (Diffusion Coefficient)
- 汤生系数 (Townsend Coefficient)
- 吸附系数 (Attachment Coefficient)

## 程序说明
程序包含三个文件：
- generate.C `gas`文件生成程序。
- read.C `gas`文件读取程序。
- CMakeLists.txt `cmake`文件，后续介绍以`cmake`编译方式为例。
- makefile `make`文件。
- ReadGas.C 通过`ROOT`程序读取`gas`文件。
- GetData.sh 通过`shell`脚本从`log`文件提取参数信息。

&emsp;&emsp;`generate.C` 程序有个 `bug` ，与 `Garfield++` 版本有关，可以在 `2.0` 版本下正常运行， `3.0` 及 `master` 版本会偶尔出现 `WARNING ENERGY OUT OF RANGE, INCREASE ELECTRON ENERGY INTEGRATION RANGE` 的中断。

&emsp;&emsp;`2.0` 版本生成的 `gas` 文件与 `3.0` 及 `master` 版本不同，不过均可通过 `read.C` 程序读取不同的物理参量绘图。

### 生成气体文件
&emsp;&emsp;气体文件的后缀为 `.gas`，创建 `result` 与 `runlog` 文件夹，`result` 用来存放 `gas` 文件，`runlog` 用来存放运行日志。运行日志包含所有的计算结果。最多可以设置 `3` 种气体，可以设置不同比分、不同气压 (atm)、不同温度 (degree centigrade)，也可修改程序设置更多气体，最多不超过 `6` 种。
```bash
mkdir build
cd build
mkdir result
make
./generate [-g1 Gas1] [-f1 Fraction1] [-g2 Gas2] [-f2 Fraction2] [-g3 Gas3] [-f3 Fraction3] [-p Pressure] [-t Temperature]
```

### 结果绘图
#### 读取 `gas` 文件
- 读取单个 `gas` 文件，利用 `read.C` 程序处理。

    ```bash
    cd build
    ./read xxx
    ```
    `xxx` 为 `gas` 文件的路径。

- 读取多个 `gas` 文件进行比较，利用 `ReadGas.C` 程序处理，可根据条件修改。对比结果保存在 `result.root` 文件中。

    ```bash
    root -l 'ReadGas.C("xxx")'
    ```
    `xxx` 为 `gas` 文件夹的路径。

#### 读取 `log` 文件
&emsp;&emsp;通过 `shell` 脚本从 `log` 文件中提取信息到 `txt` 文件，再通过 `ROOT` 程序绘图，此方法仅适用于 `3.0` 及 `master` 版本输入的 `log` 文件，`2.0` 版本输入 `log` 信息中物理参数未研究明白。


