# Garfieldpp
用于模拟计算气体探测器的物理参数，主要是GEM探测器与<sup>3</sup>He探测器。

程序编译方式官方已经全面转向`cmake`，推荐使用`cmake`，不过依旧提供`make`编译方式，需要对程序做相应调整，主要是读取的文件路径，如`ansys文件`，或者`离子数据文件`。因为程序使用`root文件`保存数据，所以`CMakeLists.txt`与`makefile`均需要链接到`ROOT`。
- 使用`make`编译，需要在环境变量`setupGarfield.sh`中添加:
    ```bash
    export LIBRARY_PATH=$GARFIELD_INSTALL/lib:$LIBRARY_PATH
    ```
- 使用`cmake`编译，需要在`CMakeLists.txt`添加:
    ```bash
    find_package(ROOT)
    ```

程序运行最好先生成所需要的气体文件(xxx.gas)，然后在程序中加载该文件，可以提高计算速度。

程序中均采用`AvalancheMicroscopic`类进行电子微观运动的描述，也可修改为`AvalancheMC`或`DriftLineRKF`类，其中`DriftLineRKF`类貌似必须加载气体文件才能正常运行。

存在的问题：
在GEM电子雪崩计算中，通过`AvalancheMC`类对离子的运动进行追踪，会存在内存泄漏的问题，暂未解决。