---
kernelspec:
  name: wolframlanguage14
  display_name: 'Mathematica'
---

# 动力系统的一些基本概念

动力系统描述一个动态演化的过程, 具体使用的数学工具是微分方程与映射, 也就是形如
$$
\dot{x}=f(x,t)
$$
与
$$
x_{n+1}=g(x_n)
$$
这样的系统, 其中 $x,x_n\in\mathbb{R}^n$, $f:\mathbb{R}^n\times \mathbb{R}\rightarrow\mathbb{R}^n$, $g:\mathbb{R}^n\rightarrow\mathbb{R}^n$ 一般默认都是具有可微性的. 当 $f,g$ 都是非线性的, 即不满足
$$
f(ax+by,t)=af(x,t)+bf(y,t),\\
g(ax+by)=ag(x)+bg(y)
$$
对任意的 $a,b\in\mathbb{R},x,y\in\mathbb{R}^n$ 成立, 那么动态演化的过程可能是非常复杂的.

我们比较熟知的力学系统一般不是 $\dot{x}=f(x,t)$ 的这种形式, 一般的是二阶微分方程的形式, 例如单自由度系统 $m\ddot{x}=f(x,t)$, 我们可以令 $y=\dot{x}$, 那么此时新的方程变为
$$
\begin{aligned}
\dot{x}&=y\\
\dot{y}&=f(x,t)/m.
\end{aligned}
$$
这样就变成了标准的一阶微分方程的形式. 对于多自由度系统的处理是类似的.

映射系统实际上也经常出现, 例如求解单变量函数的零解: $f(x)=0$, 那么最常见的方法牛顿迭代法, 实际上就是映射系统:
$$
x_{n+1}=x_n-f(x_n)/f'(x_n).
$$

```{prf:definition}
如果微分方程 $\dot{x}=f(x)$ 的左边不显含时间, 我们称这个系统是自治系统; 如果 $f$ 显含时间, 那么则称其是非自治的. 对于映射系统也有类似的概念. 如果一个迭代过程依赖于迭代的次数, 也就是可以写成 $x_{n+1}=g(x_n,n)$ 的形式, 我们称它是非自治的; 如果 $g$ 不显含迭代次数 $n$, 那么我们称其是自治的.
```


## 相空间与矢量场
对于系统 $\dot{x}=f(x,t)$ 或者 $x_{n+1}=g(x_n)$, 粗略地讲, $x$ 或 $x_n$ 所在的空间就是相空间. 通常情况下, 相空间都是 $\mathbb{R}^n$, 但也有许多情况下, 相空间可以是更复杂的空间.

```{prf:example}
考虑单摆的微分方程:
$$
\ddot{\theta}=-\frac{g}{l}\sin\theta,
$$
其中 $g$ 为重力加速度, $l$ 为摆长, $\theta\in\mathbb{T}$ 为角度. 如果我们令 $\dot{\theta}=\omega$, 那么上面的摆方程变为:

$$
\label{pendulum}
\begin{aligned}
\dot{\theta}&=\omega,\\
\dot{\omega}&=-\frac{g}{l}\sin\theta.
\end{aligned}
$$

此时的相空间为 $(\theta,\omega)$ 组成的空间, 也就是 $\mathbb{T}\times \mathbb{R}$, 它是一个柱面. 类似的, 描述双摆的微分方程, 按照这种方式写成一阶的形式, 相空间是 $\mathbb{T}^2\times \mathbb{R}^2$. 对于受约束的有限个自由度的力学系统, 相空间大多是这种曲面空间.
```

对于非自治系统 $\dot{x}=f(x,t)$, 有一种简单的方法可以将其变为自治系统. 令 $\tau=t$ 为新的时间, 那么在这个新的时间下:
$$
\begin{aligned}
\frac{dx}{d\tau}&=f(x,t),\\
\frac{dt}{d\tau}&=1,
\end{aligned}
$$
我们通过将相空间扩展一维得到了一个自治系统, 通常称这样的相空间为扩展相空间.

有了相空间的概念后, 微分方程 $\dot{x}=f(x)$ 可以认为是在每个相点 $x$ 处配备了一个切矢量 $f(x)$, 我们就把 $f$ 称作是相空间上的一个矢量场. 例如, 对于平面上的矢量场:
$$
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=-x+0.1x^3,
\end{aligned}
$$
它的矢量场的形状是这样的:

```{figure} ./asserts/figs/01矢量场.png
:alt: 图片无法加载
:width: 80%
:align: center
```

对于非自治系统, 可以认为矢量场 $f(x,t)$ 是依赖于时间的. 例如, 将刚刚的例子修改为:
$$
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=-x+0.1x^3+3\sin(2\pi t),
\end{aligned}
$$
那么矢量场随着时间变化:
```{figure} ./asserts/videos/变化矢量场.mp4
:alt: 视频无法加载
:width: 80%
:align: center
```

前面我们已经提到过, 单摆的相空间是柱面 $\mathbb{T}\times\mathbb{R}$, 那么它的矢量场如果绘制在平面上, 在水平方向上将会是周期的:
```{figure} ./asserts/figs/01摆_平面.png
:alt: 图片无法加载
:width: 80%
:align: center
```

更自然地, 我们可以将矢量场绘制在柱面上:
```{figure} ./asserts/figs/01摆_柱面.png
:alt: 图片无法加载
:width: 50%
:align: center
```
## 解与流
在通常的条件下, 给定微分方程的初值问题:
$$
\dot{x}=f(x,t),x(t_0)=x_0,
$$
解 $x(t)$ 总是存在且唯一的, 它满足上面的初值问题的条件. 为了强调解对于初值的依赖性, 很多时候将解记为 $x(t,x_0,t_0)$, 其满足:
$$
x(t_0,x_0,t_0)=x_0.
$$

如果微分方程是自治的, 那么初始的时间是无关紧要的, 即
$$
x(t,x_0,t_0)=x(t+t_1-t_0,x_0,t_1).
$$
也就是说, 只要初始位置一样, 时间间隔相同, 那么解就到达相同的位置. 这是因为, 只要 $x(t)$ 是解, 那么 $x(t+t_0)$ 也是解, 这称作是自治系统解的平移不变性.

对于自治系统, 我们一般默认初始时刻是 0 时刻, 将解记为 $x(t,x_0)=x(t,x_0,0)$.

另外, 如果微分方程的解总是存在唯一的, 那么解是可以拼接在一起的. 例如, 解在 $t_0$ 时刻在 $x_0$ 处, 在 $t_1$ 时刻在 $x_1$ 处, 那么
$$
x(t,x_0,t_0)=x(t,x_1,t_1)=x(t,x(t_1,x_0,t_1),t_1).
$$

与解相关的一个重要相关的概念是流.
```{prf:definition}
考虑微分方程 $\dot{x}=f(x,t)$, 假设其解存在唯一, 那么依赖于时间的映射 $\varphi^t$ 称作是流:
$$
\varphi^t:(x_0,t_0)\mapsto x(t,x_0,t_0).
$$
如果微分方程是自治的, 那么流 $\varphi^t$ 定义为:
$$
\varphi^t:x_0\mapsto x(t,x_0,0).
$$
```
从上面的定义可以看出, 流与解的定义十分相似, 但必须指出, 流是一个映射而不是初值问题的解. 对于自治系统, 给定时间 $t$, 流将相点映射为以该点为初值解, 经过固定时间 $t$ 以后达到的点. 从这个解释, 可以看出为什么这个映射取名为流 (英文单词为 flow). 流也被称为相流 (phase flow).

如果解存在唯一, 我们有:
$$
\varphi^s(\varphi^t(x_0,t_0),t_0+t)=\varphi^{t+s}(x_0,t_0).
$$
如果系统是自治的, 上式有个更简单的形式:
$$
\varphi^s(\varphi^t(x_0))=\varphi^{t+s}(x_0).
$$

```{prf:example}
考虑标量微分方程
$$
\dot{x}=-2x,
$$
它的通解为 $ce^{-2t}$, 那么流为 $\varphi^t(x_0)=x_0 e^{-2t}$.

考虑线性弹簧振子
$$
\label{linear_zhenzi}
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=-\omega^2x,
\end{aligned}
$$
其中 $\omega=\sqrt{k/m}$, 通解为 $c_1\cos(\omega t)+c_2\sin(\omega t)$, 那么流为:
$$
\varphi^t((x_0,y_0))=(x_0\cos(\omega t)+\frac{y_0}{\omega}\sin(\omega t),-x_0\omega\sin(\omega t)+y_0 \cos(\omega t)).
$$
```



## 解的性质

这一小节, 我们介绍常微分方程解的性质: 解的存在唯一性, 解对参数和初值的连续依赖性. 我们主要给出结论以及一些例子, 这部分偏数学的叙述, 不感兴趣的读者可跳过.

对于微分方程
$$
\dot{x}=f(x,t),
$$
在什么条件下, 给定初值 $x(t_0)=x_0$, 这个方程存在解呢? 为了严格地回答这个问题, 我们介绍一个 Lipschitz 连续的概念.
```{prf:definition}
设 $f:U\mapsto V$ 是连续映射, 其中 $U,V$ 均是 $\mathbb{R}^n$ 中的集合, 称 $f$ 是 Lipschitz 连续的, 假设存在常数 $K>0$, 使得:
$$
||f(x)-f(y)||\leq K||x-y||
$$
对任意的 $x,y\in U$ 中都成立. 假设 $U$ 是开集, $f$ 称作是局部 Lipschitz 连续的, 如果对任意的 $x\in U$, 都存在一个 $x$ 在 $U$ 中的邻域使得 $f$ 在此邻域内是 Lipschitz 连续的.
假设 $f:U\times(a,b)\mapsto \mathbb{R}^n$ 是
```

## 不变集

```{prf:definition}
对于自治系统, 相空间中的集合 $M$ 称为不变集, 若 $M$ 在相流的作用下不变, 即
$$
\forall x\in M,\forall t\in \mathbb{R}, \varphi^t(x)\in M,
$$
则称 $M$ 是不变集.
```
简单地说, 不变集就是从该集合出发的任意相点, 轨道总是在该集合内的集合. 

最简单的不变集是平衡点. 对于自治系统 $\dot{x}=f(x)$, 平衡点即满足 $f(x)=0$ 的点, 若 $x^*$ 是平衡点, 那么
$$
x(t,x^*)\equiv x^*
$$
是解. 那么单点集 $\{x^*\}$ 是一个平衡点. 对于单摆, $(\theta,\omega)=(0,0)$ 与 $(\theta,\omega)=(\pi,0)$ 都是平衡点. 对于一般的力学系统, 平衡点即系统可以在该状态下达到静平衡.

另外一个最常见的不变集是周期解. 考虑线性振子 [](#linear_zhenzi), 我们知道这个系统能量守恒:
$$
E=\frac{1}{2}m \dot{x}^2+\frac{1}{2}kx^2,
$$
也就是说 $\frac{dE}{dt}=0$. 这意味着, 在相空间中, 集合
$$
H_c=\{(x,y):m \dot{x}^2+kx^2=c\}
$$
为不变集. 当 $c=0$ 时, $H_c=\{(0,0)\}$ 为平衡点; 当 $c>0$ 时, $H_c$ 为半轴为
$$
\sqrt{\frac{c}{m}},\sqrt{\frac{c}{k}}
$$
的椭圆. 这些相空间中的椭圆实际上就是线性振子的周期解:
```{figure} ./asserts/figs/01线性振子_相图.png
:alt: 图片无法加载
:width: 80%
:align: center
```

对于保守系统, 这样通过能量积分找到不变集的方式是通用的. 对于单摆系统 [](#pendulum), 能量积分找到的周期解展现出这样的图像:
```{figure} ./asserts/figs/01摆_相图.png
:alt: 图片无法加载
:width: 80%
:align: center
```
能量积分为
$$
\{(\theta,\omega):\frac{1}{2}m l^2\omega^2-mgl\cos\theta=mgl\}
$$
对应的相轨线称为同宿轨道. 实际上, 这个能量对应于平衡点 $(\pi,0)$ 的能量, 虽然在视觉上, 连接 $(\pi,0)$ 点的相轨线是一个环路, 但是从环路上任意与 $(\pi,0)$ 不同的点出发, 要经过无穷长的时间才能达到 $(\pi,0)$, 否则, 在有限的时间达到平衡点, 则会与解的存在唯一性矛盾 (平衡点也是一个解). 这个同宿环本身也是一个不变集 (无论是否包含点 $(\pi,0)$).

除了平衡点, 周期解, 以及同宿轨道以外, 还有许多复杂的不变集. 这里我们简要介绍一下环面这个不变集. 环面相应于解就是拟周期解, 也就是类似于
$$
(\sin(\omega_1 t),\sin(\omega_2 t))
$$
的这种解. 构造这种解最简单的一个方式是拼接振子系统, 例如考虑两个线性振子的平凡拼接:
$$
\begin{aligned}
\dot{x}_1&=y_1,\\
\dot{y}_1&=-\frac{k_1}{m_1}x_1,\\
\dot{x}_2&=y_2,\\
\dot{y}_2&=-\frac{k_2}{m_2}x_2,
\end{aligned}
$$
对整体的微分方程, 两个振子的分别的能量积分为这个系统的积分. 那么单独子系统的周期解是一个一维环, 与 $\mathbb{T}^{1}$ 形状上相同, 那么在整个系统中看就是 $\mathbb{T}^{1}\times\mathbb{T}^{1}=\mathbb{T}^2$, 也就是一个环面. 这种环面的结构是通有的. 实际上, 对于这种 $n$ 自由度的 Hamilton 系统, 若系统存在 $n$ 个独立的积分, 那么系统的相空间就由 $n$ 维环面所填满 (全测度意义下, 即不是环面的集合所占体积为 0). 这个结论即 [Liouville–Arnold 定理](https://en.wikipedia.org/wiki/Liouville%E2%80%93Arnold_theorem).

除了这些规整的不变集以外, 还有混沌不变集. 著名的 Lorenz 系统:
$$
\begin{aligned}
\dot{x}& = \sigma (y - x), \\
\dot{y}& = x(\rho - z) - y, \\
\dot{z}& = xy - \beta z,
\end{aligned}
$$
当取参数为 $\sigma = 10, \quad \rho = 28, \quad \beta = \frac{8}{3}$ 时, 数值模拟可以得到蝴蝶形状的 Lorenz 吸引子
```{figure} ./asserts/figs/01Lorenz.png
:alt: 图片无法加载
:width: 60%
:align: center
```
这就是著名的混沌不变集.
## 不变集的稳定性与结构稳定性



## 分岔的粗略定义

## 离散动力系统
上述我们讨论的大多是关于微分方程的一些概念, 其中很多内容可以几乎照搬到映射系统. 例如, 对于映射 $x_{n+1}=g(x_n)$, 满足
$$
g(x)=x
$$
的点称作是不动点. 

## Poincaré 截面与扭扩
这一小节, 我们讨论微分方程与映射系统的联系. 从微分方程转变为映射的通常方法为 Poincaré 截面法.


## 混沌的含义

