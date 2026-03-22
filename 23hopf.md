# Hopf 分岔

## Hopf 分岔的基本假设
先前介绍过的鞍结分岔与音叉分岔, 分岔发生时是实轴上的特征值穿过了零. 对于二维及以上的系统, 系统平衡点处的特征值可能存在复根. 所以, 如果某对特征值在参数变化时穿越了虚轴, 那么也可能发生分岔现象.

在引言中我们简单介绍了 Hopf 分岔. 回顾其标准形式为
$$
\begin{aligned}
\dot{x}&=\mu x-y-x(x^2+y^2),\\
\dot{y}&=x+ \mu y-y(x^2+y^2).
\end{aligned}
$$
在极坐标下, 系统的半径方向的微分方程与音叉分岔的标准形式一样:
$$
\begin{aligned}
\dot{r} = &r(\mu-r^2),\\
\dot{\theta}=&1.
\end{aligned}
$$
如果我们作出 $\mu>0$ 与 $\mu<0$ 时的相图, 会得到:
```{figure} ./asserts/figs/23_hopf_phase_plot.png
:alt: 视频无法加载
:width: 100%
:align: center
```
当 $\mu<0$ 时, 原点是一个吸引的焦点; 当 $\mu>0$ 时, 原点从吸引焦点变为不稳定焦点; 同时系统具有一个**吸引的极限环** $r=\sqrt{\mu}$.
```{prf:remark}
对于自治系统, 极限环就是孤立的周期轨道. 对于无阻尼线性振子与无阻尼单摆, 周期轨道不是孤立的, 因此不是极限环. [Hillbert 第十六问题](https://en.wikipedia.org/wiki/Hilbert%27s_problems)即是, 对于平面上的多项式系统, 确定极限环的个数以及相对位置的问题. 这个问题迄今仍没有解决.
```

系统在原点处的特征值为 $\mu+i$, 可见当 $\mu$ 从负变正时, 特征值确实穿越了虚轴.

现在我们考虑更一般的形式:
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
这里我们将 $g$ 视为复域上的二元函数. $g$ 的光滑性由 $F$ 所决定. 根据多元函数的泰勒公式 (参见 [](./fulu.md)), 我们有:
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
\dot{z}=(\mu+i)z-z^2|z|.
$$
一个基本的想法是, 能否取合适的坐标变换使得方程 [](#eqz) 变成如
$$
\dot{z}=\lambda(\mu) z+\xi(\mu)z^2|z|+o(|z|^4)
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
\label{hopf_eq_norm_form}
\dot{z}=\lambda(\mu)z+c_{1}(\mu)z^2\bar{z}+o(|z|^4).
$$
系数 $c_1(0)$ 的实部为:
$$
\label{hopf_c1}
\text{Re}(c_1(0))=\frac{1}{2\omega_0}\text{Re}\Big(ig_{20}(0,0,0)g_{11}(0,0,0)+\omega_0 g_{21}(0,0,0)\Big),
$$
**其中 $g_{jk}(0,0,0)$ 为 [](#eqz) 中相应项的系数.**
```{prf:assumption}
:label: hopfa2
$\text{Re}(c_1(0))\neq 0$.
```
```{prf:theorem} Hopf 分岔定理
如果 [](#hopfa1) 和 [](#hopfa2) 成立, 当 $\text{Re}(c_1(0))< 0$ ($\text{Re}(c_1(0))> 0$), 微分方程 [](#basic_eq) 在 $\mu=0$ 附近发生超临界 (亚临界) 的 Hopf 分岔. 更具体地, 在平面上的原点附近, 对于 $\mu$ 充分小时:
- 当 $\alpha'(0)\mu<0$ 则方程的平衡点 $x=0$ 稳定; 当 $\alpha'(0)\mu>0$, 平衡点 $x=0$ 不稳定;
- 当参数满足 $\alpha'(0)\text{Re}(c_1(0))\mu<0$ 时方程存在唯一一个周期解, 周期解的振幅为 $O(|\mu|^{1/2})$; 当 $\text{Re}(c_1(0))<0$ 周期解稳定, 当 $\text{Re}(c_1(0))>0$ 周期解不稳定.
```
```{prf:proof}
先前已经论述, 在 [](#hopfa1) 的假设下, 我们将系统在原点附近进行坐标变换得到了 [](#hopf_eq_norm_form). 取极坐标 $z=re^{i\theta}$, 代入 [](#hopf_eq_norm_form) 并分别取实部与虚部得到:
$$
\dot{r}&=\alpha(\mu)r+\xi(\mu)r^3+O(|r|^5),\\
\dot{\theta}&=\omega(\mu)+\eta(\mu)r^2+O(|r|^4),
$$
其中 $c_1(\mu)=\xi(\mu)+i\eta(\mu)$. 由于在 $r=0$ 附近且 $\mu$ 较小时, 上面第二个式子的右端总是非负的, 将上下两式相除得到:
$$
\label{hopf_eq_r_theta}
\frac{dr}{d\theta}=\frac{\alpha(\mu)r+\xi(\mu)r^3+O(|r|^5)}{\omega(\mu)+\eta(\mu)r^2+O(|r|^4)}=\frac{\alpha(\mu)}{\omega(\mu)}r+\frac{\xi(\mu)\omega(\mu)-\alpha(\mu)\eta(\mu)}{\omega^2(\mu)}r^3+O(r^4),
$$
即我们得到了一个 $r$ 关于 $\theta$ 的微分方程. 假设 $\mu$ 较小且 $r$ 较小时, 方程的解由:
$$
r(\theta,r_0,\mu)=u_1(\theta,\mu)r_0+u_2(\theta,\mu)r_0^2+u_3(\theta,\mu)r_0^3+\cdots
$$
给出, 其中 $r(0,r_0,\mu)=r_0$ 为初值. 那么求导以后代回原方程 [](#hopf_eq_r_theta), 并比较同阶项, 我们得到一系列线性微分方程:
$$
\frac{du_1}{d\theta}=\frac{\alpha(\mu)}{\omega(\mu)} u_1,\quad u_1(0,\mu)=1,\\
\frac{du_2}{d\theta}=\frac{\alpha(\mu)}{\omega(\mu)} u_2,\quad u_2(0,\mu)=0,\\
\frac{du_3}{d\theta}=\frac{\alpha(\mu)}{\omega(\mu)} u_3+\frac{\xi(\mu)\omega(\mu)-\alpha(\mu)\eta(\mu)}{\omega^2(\mu)} u_1^3,\quad u_3(0,\mu)=0,
$$
连续地解上述线性微分方程 (并将得到的解代入到后续微分方程中), 我们求得:
$$
u_1(\theta,\mu)&=e^{\frac{\alpha(\mu)}{\omega(\mu)}\theta},\\
u_2(\theta,\mu)&=0,\\
u_3(\theta,\mu)&=\big(e^{\frac{3\alpha(\mu)}{\omega(\mu)}\theta}-e^{\frac{\alpha(\mu)}{\omega(\mu)}\theta}\big)\big(\frac{\xi(\mu)}{2\alpha(\mu)}-\frac{\eta(\mu)}{2\omega(\mu)}\big).
$$
因此, 有:
$$
\label{sol_r}
r(\theta,r_0,\mu)=e^{\frac{\alpha(\mu)}{\omega(\mu)}\theta}r_0+\big(e^{\frac{3\alpha(\mu)}{\omega(\mu)}\theta}-e^{\frac{\alpha(\mu)}{\omega(\mu)}\theta}\big)\big(\frac{\xi(\mu)}{2\alpha(\mu)}-\frac{\eta(\mu)}{2\omega(\mu)}\big)r_0^3+O(r_0^4).
$$

为了研究分岔出来的周期解, 我们考虑映射
$$
g(r_0,\mu)=r(2\pi,r_0,\mu),
$$
其中 $r_0>0$ 为一个小量. 如果 $g(r_0,\mu)=r_0$ 成立, 代表系统 [](#hopf_eq_norm_form) 的轨线从 $(x,y)=(r_0,0)$ 出发, 绕原点一周后再次到达 $r_0$, 即我们所寻找的周期轨道.

同研究音叉分岔时完全一样, 令
$$
H(r_0,\mu)=\frac{g(r_0,\mu)-r_0}{r_0}.
$$
由 $H$ 的定义以及 $r(\theta,r_0,\mu)$ 的表达式 [](#sol_r), $H(r_0,\mu)$ 满足:
$$
H(0,0),\quad H_{\mu}(0,0)=\frac{2\pi}{\omega_0}\alpha'(0)\neq 0.
$$
由隐函数定理, 存在函数 $\mu(r_0)$ 使得 $H(r_0,\mu(r_0))=0$, 或者等价地, $r_0=g(r_0,\mu(r_0))$. 根据这个性质, 代入 [](#sol_r) 到 $H(r_0,\mu(r_0))=0$ 可反解出:
$$
\mu(r_0)=-\frac{\omega_0}{2\pi\alpha'(0)}u_{3}(2\pi,0)r_0^2+O(r_0^3),
$$
其中
$$
u_3(2\pi,0)=\frac{2\pi\xi(0)}{\omega_0}.
$$
注意到 $\xi(0)=\text{Re}(c_1(0))$. 那么我们得到非平凡不动点:
$$
r_0(\mu)=\sqrt{-\frac{\alpha'(0)}{\xi(0)}\mu}+O(|\mu|^{3/2}).
$$

对于 $g$ 的平凡不动点 $r_0=0$, 
$$
g_{r_0}(0,\mu)=e^{\frac{\alpha(\mu)}{\omega(\mu)}}.
$$
由于 $\alpha(0)=0,\alpha'(0)\neq 0$, 所以
- 当 $\alpha'(0)\mu>0$, $g_{r_0}(0,\mu)>1$, 平凡不动点 $r_0=0$ 不稳定;
- 当 $\alpha'(0)\mu<0$, $g_{r_0}(0,\mu)<1$, 平凡不动点 $r_0=0$ 稳定.

对于 $g$ 的非平凡不动点 $r_0=r_0(\mu)$,
$$
g_{r_0}(r_0(\mu),\mu)&=e^{2\pi\frac{\alpha(\mu)}{\omega(\mu)}}+3r_0^2u_3(2\pi,\mu)+O(|r_0|^3)\\
&=1-\frac{4\pi\alpha'(0)}{\omega_0}\mu+O(|\mu|^{3/2}).
$$
由于 $\alpha'(0)\xi(0)\mu<0$ 且 $\omega_0>0$. 根据上式, 我们知:
- $\xi(0)$ 小于 $0$ 时周期解稳定;
- $\xi(0)$ 大于 $0$ 时周期解不稳定.
```
```{prf:remark}
在 $\omega_0>0$ 的假设下, $\xi(0)=\text{Re}(c_1(0))<0$ 周期解稳定, 我们称此时的分岔为超临界的 Hopf 分岔; $\xi(0)=\text{Re}(c_1(0))>0$ 周期解不稳定, 称此时的分岔为亚临界的 Hopf 分岔. 同音叉分岔一样, 参数相同时, 平凡解的稳定性总与周期解的稳定性相异.
```

## 系数 $c_{1}(\mu)$ 的计算
由上一小节的 Hopf 分岔定理, 范式的系数 $c_{1}(\mu)$ 是很关键的. 我们也给出了从最开始的复数形式的方程
$$
\dot{z}=\lambda(\mu)z+g_2(z,\bar{z},\mu)+g_3(z,\bar{z},\mu)+\cdots
$$
直接得到 $c_{1}(\mu)$ 的实部 [](#hopf_c1) 的公式. 那么, 具体应该如何计算 $c_{1}(\mu)$ 呢? 有兴趣的读者可参见 [](./61symbolic.md), 我们提供了计算的思路以及相关 Mathematica 代码.

## 投影法计算 $g_{jk}$
由 $\text{Re}(c_1(0))$ 的表达式, 由需要计算方程 [](#eqz) 中的 $g_{jk}$, 而 $g_{jk}$ 的计算涉及到对于原方程作坐标变换再求导数, 这一过程是比较繁琐的.

下面介绍的投影法可以不做坐标变换而得到 $g_{jk}$.
```{prf:observation}
如果将线性映射 $A$ 看作是复空间 $\mathbb{C}^2$ 上的线性映射, 那么它的特征向量 $q,\bar{q}$ 是 $\mathbb{C}^2$ 上的一组基, 也就是任意一个复二维向量可唯一地写为 $q,\bar{q}$ 的线性组合. 然而, 我们处理的是实向量 $x\in\mathbb{R}^2\subset\mathbb{C}^2$, 它的分解必然满足:
$$
x=zq+\bar{z}\bar{q}.
$$
如果我们想从关于 $x$ 的微分方程转变为关于 $z$ 的微分方程, 一个自然的想法是消去 $\bar{z}$ 以后求导. 将 $x$ 向与 $\bar{q}$ 垂直的矢量 $p$ 投影, 我们就得到了 $\langle x,p \rangle=z\langle q,p \rangle$. 取合适的 $p$ 使得 $\langle q,p \rangle=1$, 我们也就得到了 $z=\langle x,p \rangle$. 以上论述是下面计算的一个粗略想法.
```
如果 [](#hopfa1) 成立, 取 $q$ 为矩阵 $A$ 相应的复特征向量 (以下均省略掉对于参数 $\mu$ 的依赖性), 即 $Aq=\lambda q$. 取 $p$ 为 $A^T$ 与特征值 $\bar{\lambda}$ 相应的特征向量, 即
$$
A^{T}p=\bar{\lambda}p.
$$
取合适的 $p$ 使得 $\langle p,q\rangle=1$. 复矢量 $p$ 与 $\bar{q}$ 垂直, 即 $\langle p,\bar{q}\rangle=0$. 实际上, 
$$
\bar{\lambda}\langle p,\bar{q}\rangle&=\bar{\lambda}(\bar{p}^{T}\bar{q})=\bar{p}^{T}(\bar{\lambda}\bar{q})=\bar{p}^{T}A\bar{q},\\
\lambda\langle p,\bar{q}\rangle&=\lambda\bar{p}^{T}\bar{q}=(\lambda\bar{p}^{T}\bar{q})=\bar{p}^{T}A\bar{q}.
$$
由于 $\lambda\neq\bar{\lambda}$, 只能有 $\langle p,\bar{q}\rangle=0$.

对矢量 $x\in\mathbb{R}^2\subset\mathbb{C}^2$ 进行分解:
$$
x=zq+\bar{z}\bar{q}.
$$
两边同时与 $q$ 内积得到 $z=\langle x,q\rangle$. 对于 $z=\langle x,q\rangle$ 两边同时关于时间 $t$ 求导, 并代入 [](#basic_eq) 得到:
$$
\dot{z}&=\lambda z+\langle p,F(zq+\bar{z}\bar{q})\rangle\\
&=\lambda z+g_2(z,\bar{z})+g_3(z,\bar{z})+O(|z|^4).
$$
假设 $F$ 可以写成如下多重线性型:
$$
F(x)=\frac{1}{2}B(x,x)+\frac{1}{6}C(x,x,x)+O(||x||^4),
$$
其中
$$
B_k(\mathbf{\xi},\mathbf{\eta})=\sum_{i_1=1}^{2}\sum_{i_2=1}^{2}\frac{\partial ^2F_k}{\partial x_{i_1}\partial x_{i_2}}(0)\xi_{i_1}\eta_{i_2},\quad k=1,2\\
C_k(\mathbf{\xi},\mathbf{\eta},\mathbf{\gamma})=\sum_{i_1=1}^{2}\sum_{i_2=1}^{2}\sum_{i_3=1}^{2}\frac{\partial ^3f_k}{\partial x_{i_1}\partial x_{i_2}\partial x_{i_3}}(0)\xi_{i_1}\eta_{i_2}\gamma_{i_3},\quad k=1,2.
$$
参见 [](./fulu.md). 那么
$$
B(zq+\bar{z}\bar{q},zq+\bar{z}\bar{q})=z^2B(q,q)+2z\bar{z}B(q,\bar{q})+\bar{z}^2B(\bar{q},\bar{q}).
$$
因此
$$
g_{20}(0)&=\frac{\partial^{2}}{\partial z^2}|_{z=0,\bar{z}=0}\langle p,F(zq+\bar{z}\bar{q})\rangle\\
&=\langle p,B(q,q)\rangle.
$$
同理得到:
$$
g_{11}(0)=\langle p,B(q,\bar{q})\rangle,\\
g_{02}(0)=\langle p,B(\bar{q},\bar{q})\rangle,\\
g_{21}(0)=\langle p,C(q,q,\bar{q})\rangle
$$
## 计算实例
判断一个平面系统是否发生 Hopf 分岔, 要在分岔点验证条件 [](#hopfa1). 如果不能保证当 $\mu$ 变化时 $x=0$ 是系统的平衡解, 即平衡解依赖于参数 $\mu$, 不妨设其为 $x^{*}(\mu)$. 做坐标变换 $\xi=x-x^{*}(\mu)$, 仍然以 $\mu$ 作为参数, 即可保证新变量下 $\xi=0$ 在 $\mu$ 变动时总为平衡点.

验证条件 [](#hopfa1) 后, 则需要计算 [](#hopfa2) 中的系数 $c_{1}(0)$ 来判断其发生了 Hopf 分岔.






