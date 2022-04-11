# Garfieldpp
用于模拟计算气体探测器的物理参数，主要是GEM探测器与<sup>3</sup>He探测器。

程序编译方式官方已经全面转向cmake，推荐使用cmake，不过依旧提供make编译方式，需要对程序做相应调整，主要是读取的文件路径，如ansys文件，或者离子数据文件。

程序运行最好先生成所需要的气体文件(xxx.gas)，然后在程序中加载该文件，可以提高计算速度。

程序中均采用`AvalancheMicroscopic`类进行电子微观运动的描述，也可修改为`AvalancheMC`或`DriftLineRKF`类，其中`DriftLineRKF`类貌似必须加载气体文件才能正常运行。

存在的问题：
在GEM电子雪崩计算中，通过`AvalancheMC`类对离子的运动进行追踪，会存在内存泄漏的问题，暂未解决。