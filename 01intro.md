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

考虑单摆的微分方程:
$$
\ddot{\theta}=-\frac{g}{l}\sin\theta,
$$
其中 $g$ 为重力加速度, $l$ 为摆长, $\theta\in\mathbb{T}$ 为角度. 如果我们令 $\dot{\theta}=\omega$, 那么上面的摆方程变为:
$$
\begin{aligned}
\dot{\theta}&=\omega,\\
\dot{\omega}&=-\frac{g}{l}\sin\theta.
\end{aligned}
$$
此时的相空间为 $(\theta,\omega)$ 组成的空间, 也就是 $\mathbb{T}\times \mathbb{R}$, 它是一个柱面. 类似的, 描述双摆的微分方程, 按照这种方式写成一阶的形式, 相空间是 $\mathbb{T}^2\times \mathbb{R}^2$. 对于受约束的有限个自由度的力学系统, 相空间大多是这种曲面空间.

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
```{code-cell} wolframlanguage14
VectorPlot[{y, -x+0.1x^3}, {x, -6, 6}, {y, -6, 6}]
```
## 解与流


## 解的性质

## 不变集

## 不变集的稳定性与结构稳定性

## 分岔的粗略定义

## 离散动力系统

## Poincaré 截面与扭扩

## 混沌的含义

