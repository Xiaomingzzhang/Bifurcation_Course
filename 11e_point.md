# 平面线性系统的初等奇点分类

本章我们考虑线性系统奇点的分类, 首先我们考虑最简单的平面情况.

## 平面初等奇点的分类
考虑方程
$$
\label{linear_eq}
\dot{X}=AX,
$$
其中 $X=(x_1,x_2)^{T}$ 是二维矢量, $A$ 是 $2\times 2$ 的常数矩阵. 为了让线性系统的形式最简单, 我们通常取非奇异的变换
$$
X=BY, Y=B^{-1}X,
$$
其中 $B$ 是二阶的非奇异矩阵, $Y=(y_1,y_2)$ 也为二维矢量. 那么在新的坐标下:
$$
\dot{Y}=JY,\quad J=B^{-1}AB.
$$
如果取合适的 $B$ 使得 $J$ 为对角阵, 那么上述过程即线性代数中的对角化, 对角矩阵的对角元为 $A$ 的特征值. 然而我们知道不是所有的方阵都可以对角化, 且特征值也可能是复数. 一般地, 有下面的结论:
```{prf:theorem}
存在非奇异矩阵 $B$, 使 $J = B^{-1}AB$ 为以下形式之一 (Jordan 型):
1. $J = \begin{bmatrix} \lambda & 0 \\ 0 & \mu \end{bmatrix}$
2. $J = \begin{bmatrix} \lambda & 1 \\ 0 & \lambda \end{bmatrix}$
3. $J = \begin{bmatrix} \alpha & -\omega \\ \omega & \alpha \end{bmatrix}$
其中 $\lambda,\mu,\alpha,\omega$ 为实数, $\omega>0$.
```

```{prf:proof}
1. 如果 $A$ 的特征值都是实数, 且有两个独立的特征向量 $v_1,v_2$, 那么令矩阵 $B=[v_1,v_2]$;
1. 如果 $A$ 的二重特征值 $\lambda$ 仅有一个独立的特征向量, 那么存在 $v\neq 0$, $(A-\lambda I)v\neq 0,(A-\lambda I)^2 v=0$, 取 $B=[(A-\lambda I)v,v]$;
1. 如果 $A$ 的特征值为一对共轭复数 $\lambda=\alpha \pm i\omega$, 取 $v=v_1+i v_2$ 为 $A$ 的特征值 $\alpha - i\omega$ 对应的特征向量, 即
$$
A(v_1+i v_2)=(\alpha - i\omega)(v_1+i v_2).
$$
取 $B=[v_1,v_2]$.
```

如果 $\det A\neq0$, 则 $\lambda,\mu,\omega\neq0$, 线性方程 [](#linear_eq) 的奇点 $(0,0)$ 称为初等奇点.

以下通过 Jordan 型 1, 2 和 3 给出初等奇点的分类.

### Jordan 型为 $J = \begin{bmatrix} \lambda & 0 \\ 0 & \mu \end{bmatrix}$ 的分类

- $\lambda=\mu\neq0$, $\lambda<0$ 或 $\lambda>0$: 星形节点. 系统的流为
$$
\varphi^t(x_1,x_2)=(x_1 e^{\lambda}t,x_2 e^{\lambda t}).
$$
相轨线结构为
```{figure} ./asserts/figs/11_星形节点_不稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda>0$ 时的不稳定星形节点
```
```{figure} ./asserts/figs/11_星形节点_稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda<0$ 时的稳定星形节点
```
- $\mu<\lambda<0$ 或 $\lambda>\mu>0$: 双轴节点. 系统的流为
$$
\varphi^t(x_1,x_2)=(x_1e^{\lambda}t,x_2 e^{\mu t}).
$$
相轨线结构为
```{figure} ./asserts/figs/11_双轴节点_不稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda>\mu>0$ 时的不稳定双轴节点
```
```{figure} ./asserts/figs/11_双轴节点_稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\mu<\lambda<0$ 时的稳定双轴节点
```
- $\mu<0<\lambda$: 鞍点. 系统的流为
$$
\varphi^t(x_1,x_2)=(x_1e^{\lambda}t,x_2 e^{\mu t}).
$$
相轨线结构为
```{figure} ./asserts/figs/11_鞍点.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda>0>\mu$ 时的鞍点
```
### Jordan 型为 $J = \begin{bmatrix} \lambda & 1 \\ 0 & \lambda \end{bmatrix}$ 的分类
- 单轴节点. 系统的流为
$$
\varphi^t(x_1,x_2)=((x_1+x_2 t )e^{\lambda}t,x_2 e^{\lambda t}).
$$
相轨线结构为
```{figure} ./asserts/figs/11_单轴节点_不稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda>0$ 时的不稳定单轴节点
```
```{figure} ./asserts/figs/11_单轴节点_稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\lambda<0$ 时的稳定单轴节点
```
### Jordan 型为 $J = \begin{bmatrix} \alpha & -\omega \\ \omega & \alpha \end{bmatrix}$ 的分类
取极坐标变换 $x_1=r\cos\theta,x_2=r\sin\theta$ 得到
$$
\dot{r}=\alpha,\quad\dot{\theta}=\omega.
$$
- $\alpha\neq 0$: 焦点
```{figure} ./asserts/figs/11_焦点_不稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\alpha>0$ 时的不稳定焦点
```
```{figure} ./asserts/figs/11_焦点_稳定.png
:alt: 图片无法加载
:width: 60%
:align: center
$\alpha<0$ 时的稳定焦点
```
- $\alpha=0$: 中心
```{figure} ./asserts/figs/11_中心.png
:alt: 图片无法加载
:width: 60%
:align: center
```

假设二阶矩阵的形式为:
$$
A=\begin{pmatrix}
    a & b\\
    c & d
\end{pmatrix},
$$
$A$ 的特征值为 $\lambda_1,\lambda_2$, 那么
$$
\begin{aligned}
    p&=\text{tr}(A)=a+d=\lambda_1+\lambda_2,\\
    q&=\text{det}(A)=ad-bc=\lambda_1*\lambda_2.
\end{aligned}
$$
特征方程为:
$$
\text{det}|A-\lambda I|=\lambda^2-p\lambda+q=0
$$
因此针对 $p,q$, 容易得到如下分类:
- $q<0$, 零点为鞍点;
- $q>0$, 且 $p^2/4>q$, 零点为双轴节点;
- $q>0$, 且 $0<p^2/4=q$, 零点为单轴或星形节点;
- $q>0$, 且 $0<p^2/4<q$, 零点为焦点;
- $q>0$, 且 $p=0$, 零点为中心.
```{prf:example}
考虑有阻尼线性振子:
$$
m\ddot{x}=-kx-c\dot{x}.
$$
令 $\dot{x}=y$, 得到
$$
\begin{aligned}
    \dot{x}&=y,\\
    \dot{y}&=-\frac{k}{m}x-\frac{c}{m}y.
\end{aligned}
$$
矩阵 $A$ 相应的 $p,q$ 为:
$$
\begin{aligned}
    p&=-\frac{c}{m},\\
    q&=\frac{k}{m}.
\end{aligned}
$$
那么由于 $k,m>0,c\geq 0$, $q>0$ 总成立. 因此:
- $p=0$, 即 $c=0$ 时, 零点为中心;
- $0<p^2/4<q$, 即 $0<c<2\sqrt{km}$ 时, 零点为焦点;
- $p^2/4>q$, 即 $c>2\sqrt{km}$ 时, 零点为双轴节点;
- $0<p^2/4=q$, 即 $c=2\sqrt{km}$ 时, 零点为单轴节点.

上述四种情况分别对应于振动力学中线性振子振动的无阻尼, 小阻尼, 过阻尼与临界阻尼情况.
```

## 奇点的稳定性

以上我们仅仅是观察相图而直观得到了平面奇点的稳定性, 下面我们介绍几个稳定性的概念.

考虑自治系统:
$$
\label{sys}
\dot{x} = f(x), \quad x \in \mathbb{R}^n,
$$
其解记为 $x(t,x_0)$.

```{prf:definition} Lyapunov 稳定
对于方程 [](#sys) 的解 $x(t,x^*)$, 称其为 Lyapunov 稳定的, 如果对于任意给定的 $\epsilon > 0$, 都存在一个 $\delta = \delta(\epsilon) > 0$, 使得当 $\|x(t_0,x_0) - x(t_0,x^*)\| < \delta$ 时, 对于所有的 $t \geq t_0$, 系统的解 $x(t)$ 均满足:
$$
\|x(t,x_0) - x(t,x^*)\| < \epsilon,
$$
其中 $t_0\in\mathbb{R}$.
```

```{prf:definition} 渐进稳定
如果解 $x(t,x^*)$ 同时满足以下两个条件:
- 它是 Lyapunov 稳定的;
- 它是吸引的, 即存在 $\delta_0 > 0$, 使得当 $\|x(t_0,x_0) - x(t_0,x^*)\| < \delta_0$ 时, 有:
    $$
    \lim_{t \to \infty} \|x(t,x_0) - x(t,x^*)\| = 0,
    $$
则称解 $x(t,x^*)$ 是渐近稳定的.
```

```{prf:definition} 轨道稳定 (Poincaré 稳定)
解 $x(t,x^*)$ 称作是轨道稳定的, 如果对于任意的 $\epsilon > 0$, 都存在 $\delta > 0$, 使得当 $\|x(t_0,x_0) - x(t_0,x^*)\| < \delta_0$ 时, 对于所有 $t \geq t_0$, 都有:
$$
\rho(x(t,x_0), \mathcal{O}^{+}_{x^*}(t_0)) < \epsilon,
$$
其中距离定义为点到集合的下确界: $\rho(x, A) = \inf_{y \in A} \|x - y\|$, $\mathcal{O}^{+}_{x^*}(t_0)$ 为轨道:
$$
\mathcal{O}^{+}_{x^*}(t_0)=\{x(t,x^*):t>t_0\}.
$$
```

根据上一小节论述的线性系统的 Jordan 型, 不难证明:
- 如果平面系统 $\dot{X}=AX$ 的矩阵 $A$ 的特征值的实部全小于 $0$, 那么平衡点是稳定的, 且满足上面的三个稳定性的定义; 
- 只要存在 $A$ 的某个特征值的实部大于 $0$, 那么平衡点是不稳定的 (在上述三个稳定性的定义下);
- 如果平衡点是中心, 那么它是 Lyapunov 稳定和 Poincaré 稳定的, 而非渐进稳定.

对于平衡点, Lyapunov 稳定与 Poincaré 稳定没有区别. 但是对于周期轨道等不变集则有所区别, 例如对于单摆的周期轨道, 它们是 Poincaré 稳定而非 Lyapunov 稳定的. 这是因为由于相邻的周期轨道之间的周期是不同的, 即便有微小的偏差, 经过足够长的时间也足够产生出相位差从而使得轨线有足够的偏离.