# 课程介绍



这门课程注重分岔理论推导与计算的实践部分. 课程内容会逐步在此网站上上线. 如果存在叙述错误, 笔误以及其他意见, 请联系邮箱: xiaomingzzhang@163.com.

## 记号规定

- $\mathbb{Z}$ 表示整数集, $\mathbb{N}$ 表示非负整数集;
- $\mathbb{R}^n$ 表示 $n$ 维欧几里得空间, $x\in\mathbb{R}^n$ 表示为 $x=(x_1,\cdots,x_n)^{T}$; 特别地, $\mathbb{R}=\mathbb{R}^1$ 代表实数集;
- 不作特别说明的情况下, 向量 $x=(x_1,\cdots,x_n)$ 总默认为是一个列向量; 当需要区分列向量与行向量时, 我们将列向量写为 $(x_1,\cdots,x_n)^{T}$ 的形式;
- $\mathbb{C}^n$ 表示 $n$ 维复空间, $\mathbb{C}$ 代表复数集;
- 若 $z\in \mathbb{C}$, 或者 $q\in\mathbb{C}^n$, 那么 $\bar{z}$ 与 $\bar{q}$ 表示其共轭;
- $||\cdot||$ 表示一个矢量的欧几里得模长;
- 对于 $x=(x_1,\cdots,x_n)^{T},y=(y_1,\cdots,y_n)^{T}\in\mathbb{R}^n$, $\langle x,y \rangle$ 表示这两个矢量的内积:
  $$
\langle x,y \rangle =x^{T}y= \sum_{i=1}^{n}x_i y_i;
  $$
- 对于 $x=(x_1,\cdots,x_n)^{T},y=(y_1,\cdots,y_n)^{T}\in\mathbb{C}^n$, 内积 $\langle x,y \rangle$ 为:
  $$
\langle x,y \rangle =\bar{x}^{T}y= \sum_{i=1}^{n}\bar{x}_i y_i;
  $$
- $\mathbb{R}/a\mathbb{Z}$ 表示模 $a$ 的数集, 其中 $a$ 是正实数; 例如当 $a=2\pi$, 那么 $\mathbb{R}/2\pi\mathbb{Z}$ 就是表示模掉 $2\pi$ 的集合; 这种记号的原因是, 可以把相差 $ak$ ($k\in\mathbb{Z}$) 的数当作是同一个数; 那么 $\mathbb{R}/a\mathbb{Z}$ 也可以认为是长度为 $a$ 的圆周上的点的集合;
- 记 $\mathbb{T}=\mathbb{R}/2\pi\mathbb{Z}$ 为一维环, $\mathbb{T}^n=(\mathbb{R}/2\pi\mathbb{Z})^{n}$ 为 $n$ 维环面;
- $\{\cdot\}$ 表示集合, $\{x:\text{条件 1}\}$ 表示满足条件 1 的 $x$ 的集合.
- 记 $\mathbb{S}^n$ 为 $n$ 维空间中的球面:
$$
\{x\in\mathbb{R}^n: ||x||=1\};
$$
注意到 $\mathbb{S}^1=\mathbb{T}^1$;
- $\dot{x}$ 默认为变量 $x$ 对时间 $t$ 求导;
- 对于映射 $x_{n+1}=f(x_n)$, 也通常写为 $x'=f(x)$; $f^n(x)$ 对于映射 $f$ 而言, 表示 $f$ 的 $n$ 次复合映射作用在 $x$ 上而非 $f(x)$ 的 $n$ 次方.

## 主要参考文献

- 非线性动力学. 谢建华, 李登辉, 乐源. 科学出版社, 2018.
- Introduction to Applied Nonlinear Dynamical Systems and Chaos. S. Wiggins, Springer, 2003.
- Elements of Applied Bifurcation Theory. Yuri A. Kuznetsov. Springer, 1998.
