# GEM探测器X射线模拟
模拟5.9keV X射线在陶瓷GEM探测器中的响应，包括增益、能量分辨率、探测器原始信号与电压信号(需要知道前放的传递函数)，所有数据均保存在`result`文件夹(需提前建好)的`root`文件中，另外每个事件产生的信号会保存在`result`文件夹的`csv`文件中。

## ANSYS脚本
该脚本除了可以生成陶瓷GEM探测器的电场文件外，还生成了4个电极的权重场文件，对孔中心(Z=0)处的电场强度进行采样保存。

## 编译运行
```bash
make
./x-ray par
```
- par表示事件数

## 功能介绍
程序设置的是X射线模拟，也可以对电子、离子等粒子进行模拟计算，可以通过布尔参数设置程序功能。
- `plotField`设置是否绘制电场分布
- `plotDrift`设置是否绘制粒子漂移
- `driftIon`设置是否进行离子漂移的计算，计算存在内存泄漏问题，设置`false`
- `calSignal`设置是否进行信号计算
- `plotSignal`设置是否进行信号的绘制，包括了电子、离子、电子离子产生的信号，`calSignal`需设置为`true`