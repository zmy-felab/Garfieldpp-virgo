# 陶瓷GEM探测器增益模拟
GEM探测器的增益受多个因素影响，程序对这些因素进行模拟研究。

## 影响因素
- 陶瓷GEM膜间电压V
- 漂移区场强
- 收集区场强
- 陶瓷GEM的RIM大小
- 陶瓷GEM孔径(未添加)
- 陶瓷GEM厚度(未添加)
- 气体种类(未添加)
- 气体压强P
- 气体温度T

## 运行环境及方法
`garfieldpp` 版本为 `3.0` 或 `master` ，绘制电场线需要 `master` 版本。推荐通过 `makefile` 编译(`CMake` 编译需修改程序中 `ansys` 生成的 `lis` 文件路径以及数据保存路径及运行日志路径)，要注意`ansys`文件的放置路径与程序中设置的一致。运行结束会在`result`文件夹生成对应的`root`文件。离子的输运会导致程序内存泄漏，暂未解决。

```bash
make
./CeramicGEM -n nEvents -p Pressure -t Temperature -v Voltage -d Drift -i Induction -r Rim
```
这些参数都有缺省值，不需要全部指定。

如果需要进行不同条件下的对比，可以通过`Shell脚本`来循环运行，参考`run_induce.sh`。


## 结果分析
通过`AnaResult.C`对模拟生成的`root`文件进行遍历，对电子的运行情况进行判选。