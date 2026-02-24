# 分岔现象初探: 平衡点的单参数分岔
对于不含时间的自治系统:

$$
\dot{x}=f(x,\mu),
$$

其中 $x\in\mathbb{R}^n$ 是 $n$ 维列向量, $\mu\in\mathbb{R}$ 是一个实数, $f:\mathbb{R}^n\times\mathbb{R}\rightarrow\mathbb{R}^n$ 是一个 $n$ 维的充分光滑的矢量场. 那么系统的平衡点由以下方程来确定

$$
f(x,\mu)=0.
$$

下面我们考虑 $n=1$ 或者 $n=2$ 的情况, 来看看平衡点的基本分岔.


## 鞍结分岔

考虑如下一维方程

$$
\dot{x}=\mu+x^2=f(x,\mu).
$$

这个方程的相图随着参数的变化是这样的:
```{figure} ./asserts/videos/鞍结分岔.mp4
:alt: 视频无法加载
:width: 80%
:align: center
```
如果我们做出曲线: $\{(\mu,x):f(x,\mu)=\mu  +x^2=0\}$, 它是这样的:
```{figure} ./asserts/figs/02鞍结分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
曲线的上半支对应了当参数变化时一个不稳定的平衡点, 下半支对应一个稳定的平衡点. 并且当参数从正到负变化时, 系统在 $x=0$ 附近从无平衡点分岔出一个稳定与一个不稳定的平衡点, 这种分岔现象我们将其称为鞍结分岔.

## 跨临界分岔
```{figure} ./asserts/figs/02跨临界分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
当参数变化时, 两个平衡点的稳定性发生了交换, 这种分岔现象我们将其称为跨临界分岔.
## 音叉分岔
```{figure} ./asserts/figs/02音叉分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
从分岔图可以看出, 系统的平衡点从一个变成三个, 同时始终存在的那个平衡点的稳定性发生变化. 分岔图如同一个音叉, 因此得名音叉分岔.
## 滞后分岔
```{figure} ./asserts/figs/02滞后分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
为何这个分岔取名为滞后分岔呢? 因为当参数变化时, 系统的定性性质的改变有一个滞后效应, 到达临界点时稳定性突然变化.
## Hopf 分岔
以上我们介绍的几个分岔都是平衡点稳定性以及数目的变化. 现在我们介绍 Hopf 分岔, Hopf 分岔将从平衡点分岔出周期解. Hopf 分岔对于机翼的颤振与列车的蛇形运动是很有作用的.

考虑如下的平面微分方程:

$$
\begin{aligned}
\dot{x}&=\mu x-y-x(x^2+y^2),\\
\dot{y}&=x+ \mu y-y(x^2+y^2).
\end{aligned}
$$

这个方程在极坐标 $(r,\theta)$ 下是更简单的:

$$
\begin{aligned}
\dot{r} = &r(\mu-r^2),\\
\dot{\theta}=&1.
\end{aligned}
$$

由于在极坐标下两个方程是解耦的, 我们单独考虑第一个方程: $\dot{r} = r(\mu-r^2)=f(r,\mu)$. 

我们发现, 函数 $f(r,\mu)$ 与音叉分岔中的矢量场一致, 但此时只有 $r\geq0$ 是有意义的. 我们做出这个系统的分岔图:
```{figure} ./asserts/figs/02hopf分岔极坐标.png
:alt: 视频无法加载
:width: 80%
:align: center
```
当 $\mu>0$ 时, 此时分岔图的上半支对应着 $r=\sqrt{\mu}$ 这个平衡点.

返回到二维系统中, 这个平衡点对应于周期解:

$$
\begin{aligned}
r(t)&\equiv \sqrt{\mu},\\
\theta(t)&=\theta_0+t.
\end{aligned}
$$

如果我们还原到三维空间 $(x,y,\mu)\in\mathbb{R}^3$ 中看分岔图, 它将是这样的:
```{figure} ./asserts/figs/02hopf分岔空间.png
:alt: 视频无法加载
:width: 80%
:align: center
```