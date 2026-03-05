# 分岔现象初探

## 平衡点的单参数分岔
对于不含时间的自治系统:

$$
\dot{x}=f(x,\mu),
$$

其中 $x\in\mathbb{R}^n$ 是 $n$ 维列向量, $\mu\in\mathbb{R}$ 是一个实数, $f:\mathbb{R}^n\times\mathbb{R}\rightarrow\mathbb{R}^n$ 是一个 $n$ 维的充分光滑的矢量场. 那么系统的平衡点由以下方程来确定

$$
f(x,\mu)=0.
$$

下面我们考虑 $n=1$ 或者 $n=2$ 的情况, 来看看平衡点的基本分岔.


### 鞍结分岔

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

### 跨临界分岔
```{figure} ./asserts/figs/02跨临界分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
当参数变化时, 两个平衡点的稳定性发生了交换, 这种分岔现象我们将其称为跨临界分岔.
### 音叉分岔
```{figure} ./asserts/figs/02音叉分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
从分岔图可以看出, 系统的平衡点从一个变成三个, 同时始终存在的那个平衡点的稳定性发生变化. 分岔图如同一个音叉, 因此得名音叉分岔.
### 滞后分岔
```{figure} ./asserts/figs/02滞后分岔.png
:alt: 视频无法加载
:width: 80%
:align: center
```
为何这个分岔取名为滞后分岔呢? 因为当参数变化时, 系统的定性性质的改变有一个滞后效应, 到达临界点时稳定性突然变化.
### Hopf 分岔
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

## 平衡点的两参数分岔: 尖点分岔
以上考虑的都是单参数分岔, 双参数的分岔要复杂许多, 这里我们介绍尖点分岔.

考虑标量方程
$$
\dot{x}=c+dx-x^3=f(x,c,d).
$$
那么平衡点方程即在三维空间 $(c,d,x)\in\mathbb{R}^3$ 中的曲面
$$
\{(c,d,x): c+dx-x^3=0\}.
$$
这个曲面的形状是这样的:
```{figure} ./asserts/figs/02尖点分岔曲面.png
:alt: 视频无法加载
:width: 80%
:align: center
```
曲面中蓝色的部分表示 $\frac{\partial f}{\partial x}<0$, 红色的部分表示 $\frac{\partial f}{\partial x}>0$.

曲面的折痕处, 其法向与 $x$ 轴垂直, 即满足 $\frac{\partial f}{\partial x}=0$, 那么曲面的折痕轨线即:
$$
\{(c,d,x): c+dx-x^3=0,d-3x^2=0\}.
$$
这两个方程消去 $x$ 以后, 得到 $4d^3=27c^2$, 也就是曲面的折痕向 $(c,d)$ 平面的投影. 这个投影会得到一个尖点:
```{figure} ./asserts/figs/02尖点分岔曲线.png
:alt: 视频无法加载
:width: 80%
:align: center
```

如果我们想象一根线穿过 $(c,d)$ 平面, 那么在尖点内部的参数将有三个平衡点, 在尖点外部的参数只有一个平衡点.

音叉分岔的标准形式 $\dot{x}=dx-x^3$ 可认为具有对称性 $f(d,-x)=-f(d,x)$. 那么对这个系统施加额外的扰动:
$$
\dot{x}=dx-x^3\rightarrow\dot{x}=c+dx-x^3
$$
则破坏了这种对称性. 对称性破坏以后, 当固定 $c$ 足够小, 而让参数 $d$ 变化, 仍会出现一个平衡点分岔出三个平衡点的现象, 但此时已经没有了音叉分岔时的对称性:
```{figure} ./asserts/figs/02尖点分岔曲线_1.png
:alt: 视频无法加载
:width: 80%
:align: center
$c=1$ 时系统关于 $d$ 的分岔图
```
```{figure} ./asserts/figs/02尖点分岔曲线_2.png
:alt: 视频无法加载
:width: 80%
:align: center
$c=-1$ 时系统关于 $d$ 的分岔图
```

## 一维映射的单参数分岔
在前面的叙述中, 我们已经论述了映射系统与微分方程系统的联系. 因此, 研究映射系统的分岔同样也是有意义的. 而映射系统最简单的不变集也就是不动点, 下面我们简单考虑一维映射系统不动点的两类单参数分岔.

首先我们介绍一维区间映射系统一种方便的表示方式. 对于一维映射系统 $x_n=g(x_{n-1})$, 考虑 $g$ 的函数图像是很方便的. 首先 $g$ 的函数图像与坐标轴的对角线 $y=x$ 的交点的 $x$ 坐标即 $g$ 的不动点. 另外, 考虑一点 $x$ 在 $g$ 的映射下的轨道, 采用如下阶梯图的方式也是很方便的:
```{figure} ./asserts/figs/02阶梯图示例.png
:alt: 视频无法加载
:width: 80%
:align: center
g(x)=4x(1-x) 的阶梯图, 从 $x=0.2$ 出发, 迭代三次.
```

### 鞍结分岔



### 周期倍化分岔