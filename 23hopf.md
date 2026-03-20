# Hopf 分岔

## Hopf 分岔的基本假设
先前介绍过的鞍结分岔与音叉分岔, 分岔发生时是实轴上的特征值穿过了零. 对于二维及以上的系统, 系统平衡点处的特征值可能存在复根. 所以, 如果某对特征值在参数变化时穿越了虚轴, 那么也可能发生分岔现象.

在引言中我们简单介绍了 Hopf 分岔. 现在我们考虑更一般的形式:
$$
\label{basic_eq}
\dot{x}=f(x,\mu)=A(\mu)x+F(x,\mu),
$$
其中 $x=(x_1,x_2)\in\mathbb{R}^2,\mu\in\mathbb{R}^1$, $f$ 在 $(x,\mu)=(0,0)$ 附近充分光滑, 满足 $f(0,\mu)=0$ 和 $F(x,\mu)=o(||x||)$, 矩阵 $A(\mu)$ 的形式为

\begin{equation}
A(\mu)=f_x(0,\mu)=
\begin{pmatrix}
a(\mu)&b(\mu)\\
c(\mu)&d(\mu)
\end{pmatrix}
\end{equation}.
```{prf:assumption}
:label: hopfa1
- 当 $\mu$ 充分小时, $A(\mu)$ 有一对共轭特征值:
$$
\lambda(\mu),\bar{\lambda}(\mu)=\alpha(\mu)\pm i\omega(\mu),
$$
其中
$$
\alpha(\mu)=\frac{1}{2}\text{tr}A(\mu),\\
\omega(\mu)=\sqrt{\text{det}A(\mu)-\frac{1}{2}\alpha^2(\mu)};
$$
- $\alpha(0)=0,\omega_0=\omega(0)=\sqrt{\text{det}A(0)}>0$;
- $\alpha'(0)\neq 0$.
```
条件 $\alpha'(0)\neq 0$ 称为横截条件, 它保证当 $\mu$ 在 $0$ 附近单调变化时特征值能穿越虚轴.

首先, 我们先将系统的线性部分变为 Jordan 型. 由先前平面系统所述, 取 $\bar{\lambda}(\mu)$ 对应的特征向量 $v_1(\mu)+iv_2(\mu)$, 令 $P(\mu)=[v_1(\mu),v_2(\mu)]$. 取变换 $x=P(\mu)y$, 我们得到:
$$
\label{eqy}
\dot{y}=B(\mu)y+P^{-1}(\mu)F(P(\mu)y,\mu),
$$
其中
\begin{equation}
B(\mu)=
\begin{pmatrix}
\alpha(\mu)&-\omega(\mu)\\
\omega(\mu)&\alpha(\mu)
\end{pmatrix}.
\end{equation}
将 [](#eqy) 式的第二行乘以虚数 $i$ 再与第一式相加, 并令:
$$
z=y_1+iy_2,\quad y_1=\frac{z+\bar{z}}{2},\quad y_2=\frac{z-\bar{z}}{2i}.
$$
我们得到复数形式的方程:
$$
\label{eqz}
\dot{z}=\lambda(\mu)z+g(z,\bar{z},\mu),
$$
这里我们将 $g$ 视为复域上的二元函数. $g$ 的光滑性由 $F$ 所决定. 根据多元函数的泰勒公式, 我们有:
$$
g(z,\bar{z},\mu)=g_2(z,\bar{z},\mu)+g_3(z,\bar{z},\mu)+\cdots+g_{n}(z,\bar{z},\mu)+o(|z|^n),
$$
其中
$$
g_l(z,\bar{z},\mu)=\sum_{j+k=l}\frac{1}{j!k!}\frac{\partial^l g}{\partial z^j \partial |z|^k}(0,0,\mu)z^j\bar{z}^k, l=2,\cdots,n.
$$

## 范式方法的初步应用

对于 Hopf 分岔的标准形式:
$$
\begin{aligned}
\dot{x}&=\mu x-y-x(x^2+y^2),\\
\dot{y}&=x+ \mu y-y(x^2+y^2).
\end{aligned}
$$
如果采用刚刚的复变换, 我们得到:
$$
\dot{z}=(\mu+i)z-z|z|^2.
$$
一个基本的想法是, 能否取合适的坐标变换使得方程 [](#eqz) 变成如
$$
\dot{z}=\lambda(\mu) z+\xi(\mu)z|z|^2+o(|z|^4)
$$
的形式. 如果我们仍然取线性变换, 那么总是没有办法影响非线性项. 范式方法的主要思路就是取一系列的非线性变换使得尽可能消去系统的较低阶的项.

取坐标变换:
$$
\label{transform}
z=w+h_2(w,\bar{w},\mu),
$$
其中 $h_2$ 的形式与 $g_2$ 相同:
$$
h_2(w,\bar{w},\mu)=\sum_{j+k=l}\frac{1}{j!k!}h_{jk}(\mu)w^j\bar{w}^k,
$$
其系数 $h_{jk}$ 待定. 变换 [](#transform) 在原点附近接近于恒等变换, 它也是可逆的. 实际上, $z$ 关于 $w$ 的偏导在原点处是单位的, 根据反函数定理, 这个变换可逆.

变换 [](#transform) 的逆变换为:
$$
\label{inverse_transform}
w=z-h_2(z,\bar{z},\mu)+o(|z|^2).
$$
对 [](#inverse_transform) 关于时间 $t$ 求导, 其中 $\dot{z}$ 用 [](#eqz) 代替, 使用 [](#transform) 回到 $w$ 坐标, 即
$$
\begin{aligned}
    \dot{w}&=\dot{z}-D_z h_2(z,\bar{z},\mu)\dot{z}-D_{\bar{z}}h_2(z,\bar{z},\mu)\dot{\bar{z}}+o(|z|^3)\\
    &=\lambda z+g_2(z,\bar{z},\mu)-D_z h_2(z,\bar{z},\mu)\big(\lambda  z+g_2(z,\bar{z},\mu)\big)-D_{\bar{z}}h_2(z,\bar{z},\mu) \text{coj}\big(\lambda z+g_2(z,\bar{z},\mu)\big)+o(|z|^3)\\
    &=\lambda w +\left[\lambda h_2(w,\bar{w},\mu)-\lambda w D_z h_2(w,\bar{w},\mu)-\bar{\lambda} \bar{w} D_{\bar{z}} h_2(w,\bar{w},\mu)+g_2(w,\bar{w},\mu)\right]+o(|w|^3).
\end{aligned}
$$

二次项
$$
\left[\lambda h_2(w,\bar{w},\mu)-\lambda w D_z h_2(w,\bar{w},\mu)-\bar{\lambda} \bar{w} D_{\bar{z}} h_2(w,\bar{w},\mu)+g_2(w,\bar{w},\mu)\right]
$$
中 $w^j\bar{w}^k$ 的系数为
$$
(\lambda-j\lambda-k\bar{\lambda})h_{jk}+g_{jk},
$$
其中 $j+k=2$. 那么如果可取:
$$
h_{jk}=\frac{g_{jk}}{j\lambda-k\bar{\lambda}-\lambda}
$$
即可消去全部的二次项, 其中 $j+k=2$. 而当 $\mu=0$ 时, 所有的 $j,k$ 的组合 (使得 $j+k=2$) 都使得上式分母不为 $0$, 因此我们可消去所有的二次项.

那么在上述变换 [](#transform) 下, 新的方程为
$$
\dot{w}=\lambda w +\tilde{g}_{3}(w,\bar{w},\mu)+o(|w|^3).
$$
**注意, 消去二次项的过程会使得三次及以上的项的系数改变**. 

我们可以继续刚刚的操作, 定义合适的 $h$ 来消去三次项, 这时 $h$ 是三阶的关于 $w,\bar{w}$ 的多项式. 三阶的接近恒等映射的变换不会改变更低阶的项, 但也同样会改变四阶及以上的项, 类似的分析得到 $h_{jk}$ 的系数仍旧为:
$$
h_{jk}=\frac{\tilde{g}_{jk}}{j\lambda-k\bar{\lambda}-\lambda},\quad j+k=3.
$$
当 $\mu=0$, 在所有 $j+k=3$ 的组合中 $j=2,k=1$ 可使得分母为 $0$, 因此, 我们可消去除了 $w^2\bar{w}$ 以外的所有的三次项.

公式
$$
h_{jk}=\frac{g_{jk}}{j\lambda-k\bar{\lambda}-\lambda},\quad j+k=n,n\geq 2
$$
是通用的. 不难证明可以消去全部的四次项. 仍然使用 $z$ 的记号. 经过三次变换以后, 我们得到:
$$
\dot{z}=\lambda(\mu)z+c_{1}(\mu)z^2\bar{z}+o(|z|^4).
$$
系数 $c_1(0)$ 的实部为:
$$
\text{Re}(c_1(0))=\frac{1}{2\omega_0}\text{Re}\Big(ig_{20}(0,0,0)g_{11}(0,0,0)+\omega_0 g_{21}(0,0,0)\Big),
$$
**其中 $g_{jk}(0,0,0)$ 为 [](#eqz) 中相应项的系数.**
```{prf:assumption}
:label: hopfa2
$\text{Re}(c_1(0))\neq 0$.
```
```{prf:theorem} Hopf 分岔定理
如果 [](#hopfa1) 和 [](#hopfa2) 成立, 当 $\text{Re}(c_1(0))< 0$ ($\text{Re}(c_1(0))> 0$), 微分方程 [](#basic_eq) 在 $\mu=0$ 附近发生超临界 (亚临界) 的 Hopf 分岔. 更具体地, 在平面上的原点附近, 对于 $\mu$ 充分小时:
- 如果 $\alpha(\mu)\text{Re}(c_1(0))<0$ 则方程仅有一个平衡点 $x=0$; 当 $\alpha(\mu)\mu<0$ 平衡点 $x=0$ 稳定, 当 $\alpha(\mu)\mu>0$ 平衡点 $x=0$ 不稳定;
- 如果 $\alpha(\mu)\text{Re}(c_1(0))>0$ 则方程存在一个周期解, 周期解的振幅为 $O(|\mu|^{1/2})$; 当 $\text{Re}(c_1(0))<0$ 周期解稳定, 当 $\text{Re}(c_1(0))>0$ 周期解不稳定.
```
<!-- ## 投影法计算 $g_{jk}$
由 Hopf 分岔定理, Hopf 分岔的判定难点在于 $\text{Re}(c_1(0))$ 的计算. 由 $\text{Re}(c_1(0))$ 的表达式, 由需要计算方程 [](#eqz) 中的 $g_{jk}$, 而 $g_{jk}$ 的计算涉及到对于原方程作坐标变换再求导数, 这一过程是比较繁琐的.

下面介绍的投影法可以不做坐标变换而得到 $g_{jk}$.

## 系数 $c_{1}(\mu)$ 的计算

## 计算实例 -->




