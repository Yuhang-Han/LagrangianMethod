# 一致SGH 格式的数值分析

20260927 韩宇航

本文分析 Caramana et. al (1998) 中提出的一致交错型格式的一致性(compatible), 相容性, 守恒性. 并给出与常用总能量减动能方案的对比.



## 1. 一致性 (compatible)

首先通过与朴素离散对比, 讨论格式的一致性.

分别考虑Lagrange 描述下的Euler 方程.



### 1.1 连续性方程

考虑连续性方程
$$
\frac1{V_i}\frac{\mathrm d V_i}{\mathrm dt} = (\nabla\cdot\vec v)_i,
\tag {eq1}
$$

其中体积元 $V_i = M_i / \rho_i$, 速度$\vec v$.

Caramana 采用的方法是先用节点表示单元体积, 然后带入连续性方程反过来作为散度的逼近使用. 即文中(27) 式
$$
(\nabla\cdot \vec v )_z = \frac1{V_z} \sum_p \vec C_p^z \cdot \vec v_p,
\tag {sch1}
$$
其中$\vec C_p^z$ 为文中定义的角向量, 表示角块受到单元内压力的合力方向.

而如果采用直接离散, 则可以将方程在单元$V_z$ 上积分, 得到
$$
\frac{\mathrm d V_z}{\mathrm dt} = \int_{V_z} (\nabla \cdot \vec v) \mathrm d V
=\oint_{\partial V_z} \vec v \cdot \vec n \,\mathrm ds = \sum_i \vec v_i \cdot \vec N_i.
\tag{1.1}
$$
其中$\vec N_i = |e_i|\vec n_i$, $e_i$ 为单元的第 i 条边, $\vec n_i$ 为第i 条边上的单位法向量. $\vec v_i$ 为边的速度. 

若以节点速度逼近边的速度 $\vec v_i = (\vec v_{i1} + \vec v_{i2})$, 代入(1.1) 则可以得到 (sch1). 故Caramana 对连续性方程的离散与直接离散并无不同.



### 1.2 动量方程

对动量方程
$$
\rho \frac{\mathrm d \vec v}{\mathrm dt} = -\nabla P,
\tag {eq2}
$$
也是直接离散. 将 (2.1) 在对偶单元$V_p$ 上积分得
$$
M_p \frac{\mathrm d \vec v_p}{\mathrm dt} = -\int_{V_p} \nabla P \mathrm d V
= \oint_{\partial V_p} P \, \mathrm d \vec s = \sum_z \vec f_z^p
=: -\sum_p (\nabla_h P)_p,
\tag {sch2}
$$
其中角力$\vec f_z^p = \vec f_p^z= P_z C_p^z$.



### 1.3 内能方程

这是一致格式与朴素离散的分野.

考虑内能方程
$$
\rho \frac{\mathrm d \vec e}{\mathrm dt} = -P\nabla\cdot \vec v,
\tag {eq3}
$$
与(1.1) 使用同样的方法, 在单元上积分即可得到内能方程的朴素离散格式
$$
M_z \frac{\mathrm d e_z}{\mathrm dt} = -\frac{P_z}{V_z} \sum_p \vec C_p^z \cdot \vec v_p =: -P\sum_z (\nabla_{h'} \cdot \vec v)_z,
\tag {1.2}
$$
但这无法确保一致性, 即离散分部积分
$$
\sum_z V_zP_z(\nabla_{h'}\cdot \vec v)_z + \sum_pV_p \vec v_p\cdot (\nabla_h P)_p = 
\text {boundary work}.
\tag {1.3}
$$
继而无法保证总内能量守恒.

而Caramana 通过共同角力推导的离散格式
$$
M_z \frac{\mathrm d  e_z}{\mathrm dt} = -\sum_p \vec f_p^z\cdot  \vec v_p,
\tag {sch3}
$$
则可以做到. 详见第3节 能量守恒部分.



## 2. 相容性 (consistency)

本节的分析中暂时不考虑人工黏性和子单元压强. 

下面在二维笛卡尔坐标系 $(x,y)\in \R^2$ 下分别考虑三个方程离散形式的截断误差.

### 2.1 连续性方程

记速度分量 $\vec v = (u, v)$, 则方程(eq1) 中散度算子$\nabla\cdot \vec v = D_x u + D_y v$. 

将格式 (sch1) 中$\vec C_p^z$ 展开, 并将体积元取为单元面积$V_z = A$ 得到离散算子
$$
D_x^h u = \frac 1 {2A}\Big( (u_2-u_4)(y_2-y_4) + (u_1-u_3)(y_2-y_4)\Big)
=: \sum_{i=1}^4 w_i u_i,
\tag {2.1}
$$
其中$u_i, y_i$ 表示第i 个顶点的x 方向速度分量和y 坐标, $w_1 = (y_2-y_4)/2A$, $w_2=(y_3-y_1)/2A$, $w_3 = -w_1$, $w_4 = -w_2$. 故有性质 $\sum w_i = 0$, $w_i \sim 1/h$.

选择单元内一点$\vec x_c = (x_c, y_c)$ 作为参考点, 记$\delta x_i = x_i - x_c, \delta y_i = y_i-y_c$. 则$\sum w_i \delta y_i=0$, $\sum w_i \delta x_i=1$.

Taylor 展开
$$
u_i = u_c + \partial_x u_c \delta x_i + \partial_y u_c \delta y_i 
+ \frac12 \partial_{xx}u_c\delta x_i^2 + \partial_{xy}u_c\delta x_i\delta y_i
+ \frac12 \partial_{yy}u_c\delta y_i^2 + \mathcal O(h^3).
$$
代入(2.1) 得
$$
\begin{aligned}
D_x^h u_c &= \sum w_i u_c + 
\sum w_i\partial_x u_c \delta x_i 
+ \sum w_i\partial_y u_c \delta y_i 
\\
&+ \sum w_i\frac12 \partial_{xx}u_c\delta x_i^2 
+ \sum w_i\partial_{xy}u_c\delta x_i\delta y_i
+ \sum w_i\frac12 \partial_{yy}u_c\delta y_i^2 + \mathcal O(h^2)\\
& =\partial_x u_c + \sum w_i\frac12 \partial_{xx}u_c\delta x_i^2 
+ \sum w_i\partial_{xy}u_c\delta x_i\delta y_i
+ \sum w_i\frac12 \partial_{yy}u_c\delta y_i^2 + \mathcal O(h^2)\\
& = \partial_x u_c +\mathcal O(h).

\end{aligned}
$$
同理可以估计
$$
D_y^h u_c = \partial_y u_c + \mathcal O(h).
$$
相加得到逼近离散算子的截断误差
$$
(\nabla_h \cdot \vec v)_z = \nabla \cdot \vec v(\vec x_z) + \mathcal O (h).
\tag {2.2}
$$

### 2.2 动量方程

将角力定义代入动量方程离散格式(sch2) 得到
$$
M_p \frac {\mathrm d \vec v_p}{\mathrm d t} = \sum_z \vec f_z^p
= -\sum_z P_z \vec C_p^z.
$$
记 $\vec\xi = \vec x_z - \vec x_p$, 将单元压强展开到节点上
$$
P_z = P_p + \nabla P_p\cdot \vec \xi + \frac12 \vec \xi^T\nabla^2P_p \vec \xi + \mathcal O(h^3),
$$
代入前式得
$$
\begin {aligned}
(\nabla_hP)_p &= -\frac1{V_p}\Big(
\sum_z \vec C_p^z P_p 
+ \sum_z \vec C_p^z \nabla P_p\cdot \vec \xi 
+ \sum_z \vec C_p^z\frac12 \vec \xi^T\nabla^2P_p \vec \xi + \mathcal O(h^3)
\Big) \\
&= -\frac 1 {V_p} (\sum_z (\nabla P_p \cdot \vec \xi) \vec C_p^z) + \mathcal O(h),
\end{aligned}
$$
轮换被求和项得到
$$
\sum_z (\nabla P_p \cdot \vec \xi) \vec C_p^z = \nabla P_p \sum_z \vec\xi \otimes \vec C_p^z
$$
其中张量积应为 $-V_p I$, 但我目前还没有验证.

于是得到
$$
(\nabla_h P)_p = \nabla P_p + \mathcal O(h).
$$

### 2.3 内能方程

将角力定义代入内能方程离散格式得到
$$
M_z \frac{\mathrm d  e_z}{\mathrm dt} = -\sum_p \vec f_p^z\cdot  \vec v_p
= -P_z \sum_p \vec C_p^z \cdot \vec v_p=-P_zV_z(\nabla_h\cdot\vec v)_z,
$$
其中使用了(sch1). 左右除以 $V_z$ 得到
$$
\rho_z \frac{\mathrm d  e_z}{\mathrm dt} =-P_z(\nabla_h\cdot\vec v)_z,
$$
与原方程(eq3) 对比, 若表示$\rho_z, P_z$ 的误差不超过$\mathcal O(h)$, 则截断误差归于离散散度误差
$$
T_e = -P(\nabla_h\cdot \vec v-\nabla\cdot \vec v) = \mathcal O(h).
$$


## 3. 守恒性

本节证明格式的总动量守恒和总能量守恒.

### 3.1 总动量守恒

定义离散总动量 $\vec P_h:= \sum_p M_p \vec v_p$. 则其变化率
$$
\frac {\mathrm d \vec P_h}{\mathrm d t} = \sum_p M_p \frac {\mathrm d \vec v_p}{\mathrm d t} =\sum_p\sum_z \vec f_z^p = \sum_z \big(\sum_p \vec f_p^z \big)
$$
下面证明单元合力 $\sum_p \vec f_p^z = 0$. 将角力定义代入得到
$$
\sum_p \vec f_z^p = \sum_p P_z\vec C_z^p = P_z\sum_p \vec C_z^p = 0,
$$
因为角向量和$\vec C_z^p = 0$.

### 3.2 总能量守恒

离散总能量$E_h = \sum_z M_ze_z + \sum_p \frac12 M_p |\vec v_p|^2$. 变化率
$$
\frac {\mathrm d E_h}{\mathrm d t} = \sum_z M_z \frac {\mathrm d e_z}{\mathrm d t}
+ \sum_p M_p \frac {\mathrm d \vec v_p}{\mathrm d t}\cdot \vec v_p.
$$
由于内能方程和动量方程均通过角力离散, 将格式(sch2, sch3) 代入即可得到${\mathrm d E_h}/{\mathrm d t} =  0$.





## 4. 与其他总能量守恒方案对比

为了实现能量守恒, 更常见的做法是用守恒格式离散能量方程, 然后减去动能得到内能. 文中也讨论了一致格式的优势.

 $e = E-K$ 这种方法在动能主导时是"大数减大数", 得到的内能会有较大舍入误差.

另外, 尽管这样做可以使总能量守恒, 但是容许错误的内能-动能分配.



## 总结与未解

本文对Caramana, et al. (1998) 中提出的compatible SGH 格式进行了数值分析, 并与常见的 $e = E-K$ 方案对比.

这篇论文行文略显随性, 格式, 理论, 例子, 与评注 混合在一起, 不像我之前见过的计算数学论文那样结构化. 尽管这些数值分析论文中都有, 但是不太容易找到.

目前一个未解决的问题是, 论文中对守恒性的论述都着眼于整个计算区域上, 而没有讨论单元上的守恒. 使用角块功逼近内能方程, 单元上的能量守恒应该可以做到. 但是证明总动量守恒中交换了求和号, 这在单元上恐怕做不到.
