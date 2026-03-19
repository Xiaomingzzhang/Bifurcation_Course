# 鞍结分岔

在本节我们叙述一个一般的一维自治系统在什么条件下会发生音叉分岔. 在叙述之前我们先来回顾一下隐函数定理.

## 隐函数定理回顾

```{prf:theorem}隐函数定理
设 $f : \mathbb{R}^{n+m} \to \mathbb{R}^m$ 为一个连续可微函数, 并设 $\mathbb{R}^{n+m}$ 具有坐标 $(\mathbf{x}, \mathbf{y})$. 固定一点 $(\mathbf{a}, \mathbf{b}) = (a_1, \dots, a_n, b_1, \dots, b_m)$ 使得 $f(\mathbf{a}, \mathbf{b}) = \mathbf{0}$, 其中 $\mathbf{0} \in \mathbb{R}^m$ 是零向量. 

如果雅可比矩阵:

$$J_{f, \mathbf{y}}(\mathbf{a}, \mathbf{b}) = \left[ \frac{\partial f_i}{\partial y_j}(\mathbf{a}, \mathbf{b}) \right]$$

是可逆的, 那么存在一个包含 $\mathbf{a}$ 的开集 $U \subset \mathbb{R}^n$, 使得存在唯一的函数 $g : U \to \mathbb{R}^m$, 满足 $g(\mathbf{a}) = \mathbf{b}$, 且对于所有 $\mathbf{x} \in U$, 都有 $f(\mathbf{x}, g(\mathbf{x})) = \mathbf{0}$. 

此外, $g$ 是连续可微的. 将前一节所示雅可比矩阵的左半部分记为:

$$J_{f, \mathbf{x}}(\mathbf{a}, \mathbf{b}) = \left[ \frac{\partial f_i}{\partial x_j}(\mathbf{a}, \mathbf{b}) \right]$$

则 $g$ 在 $U$ 中的偏导数雅可比矩阵由以下矩阵乘积给出:

$$\left[ \frac{\partial g_i}{\partial x_j}(\mathbf{x}) \right]_{m \times n} = -[J_{f, \mathbf{y}}(\mathbf{x}, g(\mathbf{x}))]_{m \times m}^{-1} [J_{f, \mathbf{x}}(\mathbf{x}, g(\mathbf{x}))]_{m \times n}.$$
```

当 $n=m=1$, 上面的结论就退化为: 如果 $f(x,y)=0$ 是平面上的光滑曲线, 且在某个点 $d_y f(x,y)\neq 0$, 那么这个在这个点附近, 曲线可以表示为一个函数 $y=g(x)$ 的图像, 即 $f(x,g(x))=0$ 在这个点附近局部地成立.

例如球面 $\mathbb{S}^2=\{(x,y,z): ||(x,y,z)||=1\}$ 可以写成是函数 $f(x,y,z)-1=0$, 其中 $f(x,y,z)=||(x,y,z)||$. 当 $z>0$ 时, $d_z f(x,y,z)=2z>0$, 那么由隐函数定理, 球面可以写成是函数图像 $z=g(x,y)$, 实际上, 我们知道 $g(x,y)=\sqrt{1-x^2-y^2}$.

## 鞍结分岔定理

```{prf:assumption}
:label: a1
假设一维矢量场 $\dot{x}=f(x,\mu)$ 是 $C^2$ 类的双变量函数, 且
$$
\label{hyp_bif}
f(0,0)=0, \quad f_{x}(0,0)=0;\\
f_\mu(0,0)=a\neq 0, \quad f_{xx}(0,0)=2b\neq 0.
$$
```
首先注意到, [](#hyp_bif) 中的前两个等式实际上是要求平衡点发生分岔, 后面两个条件才是刻画鞍结分岔发生的关键所在.

如果 [](#a1) 满足, 那么 $f$ 在 $(0,0)$ 附近可以展开为:
$$
f(x,\mu)=a\mu+bx^2+o((|\mu|,|x|^2)).
$$
首先我们来看截断矢量场 $g(x,\mu)=a\mu+bx^2$ 的分岔图:
```{figure} ./asserts/figs/21_鞍结.png
:alt: 视频无法加载
:width: 100%
:align: center
```
我们可以看到, 尽管 $a,b$ 正负不同可能导致抛物线的开口方向以及上下两支平衡点稳定性的交换, 总体上当参数变化时系统仍然是分岔出两个稳定性相异的平衡点.

由于 $f_{\mu}(0,0)\neq 0$, 那么根据隐函数定理, 存在一个小的开区间 $I$, 使得我们在 $I$ 上可以反解出 $\mu(x)$ 使得
$$
f(x,\mu(x))=0.
$$
利用 $f$ 的展开式, 可解出
$$
\mu(x)=-\frac{b}{a}x^2+o(x^2).
$$
上述函数的图像即在 $(0,0)$ 附近的平衡点曲线. 上述表达式意味着, 在 $(\mu,x)$ 平面上, 平衡点曲线与抛物线 $\mu(x)=-\frac{b}{a}x^2$ 在 $(0,0)$ 点是二次相切的. 因此, 同抛物线 $\mu(x)=-\frac{b}{a}x^2$ 一样, 方程 $f(x,\mu)$ 当 $ab\mu>0$ 时无平衡点, 当 $\mu=0$ 时一个平衡点, $ab\mu<0$ 时两个平衡点. 这两个平衡点为:
$$
x_{\pm}(\mu)=\pm\sqrt{-a\mu/b}+o(|\mu|^{1/2}).
$$

接下来我们看当 $ab\mu<0$ 时两个平衡点的稳定性, 我们有
$$
\frac{\partial f}{\partial x}(x_{\pm}(\mu),\mu)=2bx_{\pm}(\mu)+o(|\mu|^{1/2}).
$$
那么有以下情况:
- $b>0$ 时 $x_{-}(\mu)$ 稳定, $x_{+}(\mu)$ 不稳定;
- $b<0$ 时 $x_{-}(\mu)$ 不稳定, $x_{+}(\mu)$ 稳定.
  
上述定性分类与截断方程 $\dot{x}=a\mu+b\mu^2$ 一致. 那么我们有下面的定理
```{prf:theorem}鞍结分岔定理
假设矢量场 $f(x,\mu)$ 满足 [](#a1). 那么方程 $\dot{x}=f(x,\mu)$ 发生鞍结分岔. 更具体地, 对于充分小的 $\mu$, 
- 当 $ab\mu>0$ 微分方程没有平衡点;
- 当 $ab\mu<0$ 微分方程有两个稳定性互异的平衡点 $x_{\pm}(\mu)$, 其距离原点的距离为 $O(|\mu|^{1/2})$.
```
## 平面矢量场的鞍结分岔: 降阶方法初探

我们通过一个例子来考察这个问题. 考虑有阻尼的常力矩作用下的摆方程:
$$
\ddot{\theta}+\dot{\theta}+\sin\theta =M.
$$
令 $\theta=x,\dot{\theta}=y$, 上面的方程写成二阶形式:
$$
\label{eqpendulum}
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=-y-\sin(x)+M.
\end{aligned}
$$
令上面方程的矢量场为零得到 $M-\sin(x)=0$. 显然在 $M=1$ 附近, 当 $M>1$ 时无平衡点, 当 $M<1$ 时有两个平衡点. 因此 $M=1$ 时系统发生了分岔现象.
```{figure} ./asserts/figs/21_摆_鞍结.png
:alt: 视频无法加载
:width: 100%
:align: center
```
为了转化为熟悉的形式, 我们对参数与变量进行平移处理使得分岔发生在参数与状态坐标的零点. 令
$$
M=1+\mu,\quad x=z_1+\frac{\pi}{2},\quad y=z_2.
$$
那么在新变量下, [](#eqpendulum) 变为:
$$
\begin{aligned}
\dot{z}_1&=z_2,\\
\dot{z}_2&=-z_2-\cos(z_1)+1+\mu=-z_2+\mu+\frac{z_1^2}{2}+o(z_1^4).
\end{aligned}
$$
当 $\mu=0$ 时, 线性部分为:
\begin{pmatrix}
    0&1\\
    0&-1
\end{pmatrix}
取坐标变换:
\begin{equation}
\begin{pmatrix}
    z_1\\
    z_2
\end{pmatrix}=
\begin{pmatrix}
    1&-1\\
    0&1
\end{pmatrix}
\begin{pmatrix}
    x_1\\
    x_2
\end{pmatrix}.
\end{equation}
那么新变量下的方程为:
$$
\begin{aligned}
\dot{x}_1&=-\cos(x_1-x_2)+1+\mu,\\
\dot{x}_2&=-x_2-\cos(x_1-x_2)+1+\mu.
\end{aligned}
$$
由于我们是在平衡点 $(0,0)$ 以及 $\mu=0$ 附近考虑问题, 将上面的方程在 $(0,0)$ 附近展开得到
$$
\label{final_eq}
\begin{aligned}
\dot{x}_1&=\mu+\frac{1}{2}(x_1-x_2)^2+o(||(x_1,x_2)||^2),\\
\dot{x}_2&=-x_2+\mu+\frac{1}{2}(x_1-x_2)^2+o(||(x_1,x_2)||^2).
\end{aligned}
$$

现在我们实施降阶. 注意到方程 [](#final_eq) 的第二个式子的线性部分是非退化的. 那么考虑:
$$
g(x_1,x_2,\mu)=-x_2+\mu+\frac{1}{2}(x_1-x_2)^2+o(||(x_1,x_2)||^2)=0.
$$
由 $g_{x_2}(0,0,0)=-1$, 根据隐函数定理, 我们可得到:
$$
x_2(x_1,\mu)=-\mu+o(x_1^2+\mu^2).
$$
将 $x_2(x_1,\mu)$ 代入到第一个方程中, 我们得到:
$$
\dot{x}_1&=\mu+\frac{1}{2}(x_1+\mu)^2+o(x_1^2+\mu^2)=f(x_1,\mu).
$$
这是一维方程. 那么我们可以验证假设 [](#hyp_bif) 的条件:
$$
f(0,0)=0,\quad f_{x_1}(0,0)=0,\\
f_\mu(0,0)=1,\quad f_{x_{1}x_{1}}=1.
$$
所以我们断定原系统在 $(\pi/2,0)$ 以及 $M=1$ 处发生了鞍结分岔.