本次代码为matlab代码

main.m -- 主文件运行代码，其中主要参数含义为如下:
    

   
      alpha -- Netwon迭代中的alpha
        pde -- 包含真解u和右端项f. 可选example1, example2, example3. 分别对应三个例子
        g_D -- 边界值
          f -- 右端项f
      delta -- 正则化Netwon法中的delta参数
  delta_max -- 逐渐减小delta的正则化Netwon法中初始的delta参数
     method -- 迭代方法. 可选'Euler'，'Netwon'，'Adaptive_Newtow_delta'，分别对应显示迭代法，              
               Netwon法，逐渐减小delta的正则化Netwon法.
      Ftype -- F[u]中F的类型.可选'MAWS'，'MAWS_delta'，分别对应非正则化和正则化类型.
        说明 -- 当method选'Adaptive_Newtow_delta'时, delta在后面重新初始赋值delta_max，Ftype在后   
               面重新初始赋值'MAWS_delta'
         Nn -- 向量，后面对N赋值用。N=Nn(in),表示本次中网格剖分尺寸 h=1/N
     WideNn -- 向量，后面对WideN赋值用。WideN=WideNn(iwiden),表示本次中G_theta取模版宽度为WideN

     Err，K，Time，Delta -- 统计参数.分别表示对应h = 1/N, 模版宽度取WideN时数值解与真解的Inf误差，迭代
                           次数，计算所用时间，以及用'Adaptive_Newtow_delta'时，最终delta的大小.



Delta_e.m -- 报告中算法2的实现，输出Delta_eu, Gradu_Deltaeu分别对应u向方向e求二阶差分和对应          
             的梯度u(数组表示).

Eulerexplicit.m -- 显示迭代法的程序.

example1.m -- 报告中第一个例子的程序.

example2.m -- 报告中第二个例子的程序.

example3.m -- 报告中第三个例子的程序.

G_theta.m -- 报告中算法1中的实现.

GlobaltoLocalidx.m -- 拉升后的列向量对应位置s到二维数组中的[i,j]值.

Grad_MAWS_h_theta.m -- 非正则化Netwon法的梯度矩阵组装.

Grad_MAWS_h_theta_delta.m -- 正则化Netwon法的梯度矩阵组装.

Grad_uOfDelta_eu.m -- 报告中算法3的实现，对Delta_eu关于梯度u的矩阵实现.(稀疏矩阵)

initial.m -- 初始赋值为x^2+y^2的函数.

LocaltoGlobalidx.m -- 二维数组中的[i,j]拉升为列向量后对于的s位置.

MAWS_h_theta.m -- 非正则化版本的F.

MAWS_h_theta_delta.m -- 正则化版本的F.

max_delta.m -- max^{delta}函数.

min_delta.m -- min^{delta}函数.

Newton_method.m -- Netwon法程序.

Partial_xmax_delta.m -- max^{delta}函数关于x方向求偏导.

Partial_xmin_delta.m -- min^{delta}函数关于x方向求偏导.

U_plus.m -- max{u,0}函数.

U_plus_delta.m -- max^{delta}{u,0}函数.