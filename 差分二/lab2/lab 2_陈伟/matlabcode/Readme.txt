   本次代码为matlab代码,点击main.m文件，设置mian.m中Initialize的主要参数，再运行main.m即可.各文件以及其中参数含义如下. 程序会调用ifem!!!

main.m -- 主文件运行代码，其中Initialize中主要参数含义为如下:
    
          N -- 网格剖分数
          n -- 当前网格剖分数，也有h=2/n
        pde -- 包含真解u和右端项f. 可选example1, example2, example3. 分别对应三个例子
       node -- 所有网格点的横纵坐标
         N0 -- 当前计算的内点数
         u0 -- 当前计算的数值解
          P -- Wp2范数中的p
      R,Lam -- 赋初值中的R与Lambda
        Tol -- 终止条件

Conhull.m -- 计算当前点I梯度凸包面积.

Delta_e.m -- Wp2范数中的Delta

example1.m -- 报告中第一个例子的程序.

example2.m -- 报告中第二个例子的程序.

example3.m -- 报告中第三个例子的程序.

Find_delta.m -- 找使得次微分凸包面积等于fi的值

Flip.m -- 边的flip函数

Inducled_mesh.m -- 提升边界值后要Inducle mesh的函数

isconve.m -- 判断当前边连接的两个网格是否是凸的

P2.m -- 初始赋值的二次函数

Thresold.m -- 找当前点阈值的函数

Wp2err.m -- Wp2误差的计算函数

Print.m -- 打印最后结果的函数

Jump.m -- 计算边跳量的函数

example1result.mat, example2result.mat, example3result.mat -- 三个example算例的结果
