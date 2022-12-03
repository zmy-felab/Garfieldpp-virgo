# GEM电子透过率计算

GEM的电子透过率由`入孔系数`与`出孔系数`决定，是一个十分重要的参数。程序采用`Garfield++`提供的三种电荷输运的类(只能选其一)中的`DriftElectron()`函数对电子在GEM探测器中的透过率进行模拟，电子的起始与结束位置保存在`root`文件中。

程序中布尔参数的设置：
- `useRKF`表示是否使用`DriftLineRKF`类
- `useMc`表示是否使用`AvalancheMC`类
- `useMicro`表示是否使用`AvalancheMicroscopic`类

三种电荷输运类均可以获得漂移线上包含的点，可以追踪每一点的信息。(因为结果不需要，程序未提供该功能)
