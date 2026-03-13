# $n$ 阶线性系统

这一小节我们利用矩阵的 Jordan 型理论, 给出线性自治系统通解的结构. 并介绍判断平衡点稳定的 Routh-Hurwitz 准则. 最后, 给出在非线性项的扰动下, 在平衡点附近保持结构稳定性的条件.

## 解的结构
考虑 $n$ 维线性系统
$$
\label{linear}
\dot{x}=A\dot{x},
$$
其中 $A$ 为 $n\times n$ 阶常数矩阵, $x\in\mathbb{R}^n$ 为 $n$ 维矢量.
定义矩阵指数:
$$
e^{A} =\sum _{k=0}^{\infty }A^{k} /k!.
$$
上面的级数对于矩阵范数来说是收敛的, 主要原因是 $\|A^{k}\|\leq \|A\|^k$.

那么方程 [](#linear) 的通解为 $x(t)=e^{At}x_0$. 实际上, 可以形式地验证:
$$
\dot{x}=\frac{d}{dt} (I+tA+\cdots +\frac{t^{k} }{k!} A^{k} +\cdots )x_{0} =A(e^{At} x_{0} ).
$$

设 $A$ 有实特征值 $\lambda _{j} (j=1,\cdots ,k)$ 和复特征值 $\lambda _{j} =\alpha _{j} +i\omega _{j} ,\bar{\lambda }_{j} =\alpha _{j} -i\omega _{j} (j=k+1,\cdots ,n)$. 由线性代数理论可知, 存在与实特征值 $\lambda _{j}$ 对应的广义特征向量 $u_{j} (j=1,\cdots ,k)$ 和与 $\bar{\lambda } _{j} =\alpha _{j} -i\omega _{j}$ 对应的广义特征向量 $w_{j} =u_{j} +iv_{j} (j=k+1,\cdots ,n)$, 这些广义特征向量构成了 $R^{n}$ 的一组基, 并且由这组基构成的非奇异矩阵
$$
P=\left[u_{1} ,\cdots ,u_{k} ,v_{k+1} ,u_{k+1} ,\cdots v_{n} ,u_{n} \right],
$$
使得
$$
J=P^{-1} AP=\left[\begin{array}{ccc} {B_{1} } & {} & {} \\ {} & {\ddots } & {} \\ {} & {} & {B_{r} } \end{array}\right], 
$$

其中 $B=B_{j} (j=1,\cdots,r)$ 是实的方阵. 如果 $\lambda$ 是 $A$ 的实特征值之一, 则
$$
\label{type1}
B=\left[\begin{array}{cccc} {\lambda } & {1} & {0} & {0} \\ {0} & {\lambda } & {1} & {0} \\ {0} & {0} & {\ddots } & {1} \\ {0} & {0} & {0} & {\lambda } \end{array}\right], 
$$
如果 $\lambda =\alpha +i\omega$ 是 $A$ 的复特征值之一, 则
$$
\label{type2}
B=\left[\begin{array}{ccccc} {D} & {I_{2} } & {O} & {\cdots } & {O} \\ {O} & {D} & {I_{2} } & {\cdots } & {O} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } \\ {O} & {\cdots } & {O} & {D} & {I_{2} } \\ {O} & {\cdots } & {O} & {O} & {D} \end{array}\right], 
$$
其中
$$
D=\left[\begin{array}{cc} {\alpha } & {-\omega } \\ {\omega } & {\alpha } \end{array}\right],I_{2} =\left[\begin{array}{cc} {1} & {0} \\ {0} & {1} \end{array}\right],O=\left[\begin{array}{cc} {0} & {0} \\ {0} & {0} \end{array}\right]
$$

对方程 [](#linear) 进行非奇异的坐标变换 $x=Py$, 则得到 $\dot{y}=Jy$, 其中 $J=P^{-1}AP$.
在新变量下, 解为 $y(t)=e^{Jt} y_{0}$. 由矩阵指数的定义, 我们有
$$
x(t)=Pe^{Jt} P^{-1} x_{0}. 
$$ 

以下分析矩阵
$$
e^{Jt} =\left[\begin{array}{ccc} {e^{B_{1} t} } & {\cdots } & {O} \\ {O} & {\ddots } & {O} \\ {O} & {O} & {e^{B_{r} t} } \end{array}\right]
$$
的结构．如果 $B_{j} =B$ 是 $m\times m$ 阶的 [](#type1) 形式的矩阵, 则
$$
e^{Bt} =e^{\lambda t} e^{Nt} =e^{\lambda t} \left[\begin{array}{ccccc} {1} & {t} & {t^{2} /2!} & {\cdots } & {t^{m-1} /(m-1)!} \\ {0} & {1} & {t} & {\cdots } & {t^{m-2} /(m-2)!} \\ {0} & {0} & {1} & {\cdots } & {t^{m-3} /(m-3)!} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {t} \\ {0} & {0} & {0} & {\cdots } & {1} \end{array}\right], 
$$
其中
$$
N=\left[\begin{array}{ccccc} {0} & {1} & {0} & {\cdots } & {0} \\ {0} & {0} & {1} & {\cdots } & {0} \\ {0} & {0} & {0} & {\cdots } & {0} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {1} \\ {0} & {0} & {0} & {\cdots } & {0} \end{array}\right]
$$
是 $m$ 阶零幂矩阵, 即满足
$$
N^{2} =\left[\begin{array}{cccccc} {0} & {0} & {1} & {0} & {\cdots } & {0} \\ {0} & {0} & {0} & {1} & {\cdots } & {0} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } & {1} \\ {0} & {0} & {0} & {0} & {\cdots } & {0} \\ {0} & {0} & {0} & {0} & {\cdots } & {0} \end{array}\right],\cdots ,
N^{m-1} =\left[\begin{array}{cccccc} {0} & {0} & {0} & {0} & {\cdots } & {1} \\ {0} & {0} & {0} & {0} & {\cdots } & {0} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } & {0} \\ {0} & {0} & {0} & {0} & {\cdots } & {0} \\ {0} & {0} & {0} & {0} & {\cdots } & {0} \end{array}\right],
N^{k} =O,k\ge m.
$$  

现设 $B_{j} =B$ 是 $2m\times 2m$ 阶 [](#type2) 形式的矩阵．首先, 对 [](#type2) 中的 $2m\times 2m$ 矩阵 $D$,
$$
e^{Dt} =e^{\alpha t} \left[\begin{array}{cc} {\cos \omega t} & {-\sin \omega t} \\ {\sin \omega t} & {\cos \omega t} \end{array}\right]=e^{\alpha t} R, 
$$
其中
$$
R=\left[\begin{array}{cc} {\cos \omega t} & {-\sin \omega t} \\ {\sin \omega t} & {\cos \omega t} \end{array}\right], 
$$
于是
$$
e^{Bt} =e^{\alpha t} diag\left\{R\right\}e^{Nt} =e^{\alpha t} \left[\begin{array}{ccccc} {R} & {Rt} & {Rt^{2} /2} & {\cdots } & {Rt^{m-1} /(m-1)!} \\ {O} & {R} & {Rt} & {\cdots } & {Rt^{m-2} /(m-2)!} \\ {O} & {O} & {R} & {\cdots } & {Rt^{m-3} /(m-3)!} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {Rt} \\ {O} & {O} & {O} & {\cdots } & {R} \end{array}\right], 
$$
其中
$$                 
N=\left[\begin{array}{ccccc} {O} & {I_{2} } & {O} & {\cdots } & {O} \\ {0} & {O} & {I_{2} } & {\cdots } & {O} \\ {O} & {O} & {O} & {\cdots } & {O} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {I_{2} } \\ {O} & {O} & {O} & {\cdots } & {O} \end{array}\right],
$$
是 $2m$ 阶零幂矩阵, 满足

$$
N^{2} =\left[\begin{array}{cccccc} {O} & {O} & {I_{2} } & {O} & {\cdots } & {O} \\ {O} & {O} & {O} & {I_{2} } & {\cdots } & {O} \\ {O} & {O} & {O} & {O} & {\cdots } & {O} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } & {I_2} \\ {O} & {O} & {O} & {O} & {O} & {O} \end{array}\right],N^{m-1} =\left[\begin{array}{cccccc} {O} & {O} & {O} & {O} & {\cdots } & {I_{2} } \\ {O} & {O} & {O} & {O} & {\cdots } & {O} \\ {O} & {O} & {O} & {O} & {\cdots } & {O} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } & {O} \\ {O} & {O} & {O} & {O} & {O} & {O} \end{array}\right],N^{k} =O,k\ge m.
$$

根据方程 [](#linear) 的通解结构, 可得如下定理:
```{prf:theorem}
线性系统 [](#linear) 的解是形如 $t^{k} e^{\alpha t} \cos \omega t$ 和 $t^{k} e^{\alpha t} \sin \omega t$ 函数的线性组合, 其中 $\alpha$ 或 $\omega$ 之一可为零, 或两者同时为零, $k$ 为非负整数, $k+1$ 的最大值等于 $A$ 的 Jordan 型中最大 Jordan 块的阶数.
```

## 不变子空间与稳定性

考虑系统 [](#linear),设 $A$ 的特征值 $\lambda _{j} =\alpha _{j} +i\omega _{j}$ 对应的广义特征向量为 $w_{j} =u_{j} +iv_{j}$ (如果 $\omega _{j} =0$, 则相应有$v_{j} =0$)．令
$$
B=\left\{u_{1} ,\cdots ,u_{k} ,u_{k+1} ,v_{k+1} ,\cdots ,u_{m} ,v_{m} \right\},
$$
是 $R^{n}$ 的一组基 ($n=2m-k$)．定义
$$
\begin{array}{c} {E^{s} =Span\left\{u_{j} ,v_{j} \left|\alpha _{j} <0\right. \right\},} \\ {E^{c} =Span\left\{u_{j} ,v_{j} \left|\alpha _{j} =0\right. \right\},} \\ {E^{u} =Span\left\{u_{j} ,v_{j} \left|\alpha _{j} >0\right. \right\},} \end{array} 
$$
即 $E^{s}, E^{c}$ 和 $E^{u}$ 分别是实部小于 0, 等于 0 和大于 0 的特征值对应的广义特征向量 (实部和虚部) 张成的线性子空间．$E^{s},E^{c}$ 和 $E^{u}$ 分别称为系统 [](#linear) 的稳定, 中心和不稳定线性子空间.

```{prf:theorem}
在相流的作用下, $E^{s},E^{c}$ 和 $E^{u}$ 都是不变的子空间.
```
```{prf:proof}
仅考虑 $E^{s}$ 的情况. 首先证明 $E^s$ 在 $A$ 的作用下不变. 设 $A$ 的特征值为 $\lambda_1, \dots, \lambda_k$, 其中 $\text{Re}(\lambda_i) < 0$ 的特征值对应的代数重数为 $n_i$. 稳定子空间定义为这些特征值对应的广义特征子空间的直和:
$$
E^s = \bigoplus_{\text{Re}(\lambda_i) < 0} \text{Ker}(A - \lambda_i I)^{n_i}.
$$

只需证明对于 $E^s$ 中的任意一个广义特征子空间 $V_i = \text{Ker}(A - \lambda_i I)^{n_i}$, 若 $v \in V_i$, 则 $Av \in V_i$.

若 $v \in V_i$, 则满足 $(A - \lambda_i I)^{n_i} v = 0$. 而我们有:
$$
(A - \lambda_i I)^{n_i} (Av) = A \left[ (A - \lambda_i I)^{n_i} v \right] = A(0) = 0,
$$
其中我们用到了算子 $A$ 与自身生成的任何多项式算子 (如 $(A - \lambda_i I)^{n_i}$) 都是可交换的.
这说明 $Av \in \text{Ker}(A - \lambda_i I)^{n_i}$, 即 $Av \in V_i \subseteq E^s$.

如果 $E^s$ 在 $A$ 的作用下不变, 那么通过数学归纳法易证 $E^s$ 在 $A^k$ 的作用下对于任意 $k \in \mathbb{N}$ 也是不变的, 即若 $v \in E^s$, 则 $A^k v \in E^s$.
对于有限项的和 $S_N v = \sum_{k=0}^{N} \frac{t^k}{k!} A^k v$, 由于 $E^s$ 是线性子空间 (对加法封闭), 每一项 $\frac{t^k}{k!} A^k v$ 都在 $E^s$ 中, 因此其有限和 $S_N v \in E^s$.
当 $N \to \infty$ 时, 级数收敛到 $e^{At} v$. 因为子空间是闭的, 极限点也必然属于该空间:
   $$
   e^{At} v = \lim_{N \to \infty} S_N v \in E^s.
   $$
```
```{figure} ./asserts/figs/12_子空间.png
:alt: 视频无法加载
:width: 80%
:align: center
```

```{prf:theorem}
- 如果矩阵 $A$ 的所有特征值具有负实部, 那么方程 [](#linear) 的零解 $x=0$ 是渐近稳定的;
- 如果矩阵 $A$ 的某些特征值具有零实部, 但是与这些特征值对应的初等因子均是一阶的, 而其它的特征值具有负实部, 那么 $x=0$ 是稳定的, 但非渐近稳定的;
- 除上述两种情况外, 平衡点 $x=0$ 是不稳定的.
```

上述定理给出了线性系统稳定性的完整判据. 在科学研究中我们自然最关注稳定的情况, 即矩阵 $A$ 的特征值全都有负实部. 下面的 Routh-Hurwitz 准则给出了 $A$ 的特征值全都有负实部的充分必要条件.

方程 [](#linear) 零解为渐近稳定性的充要条件是特征方程
$$
\label{chara}
\det (\lambda I-A)=0 
$$
的所有根具有负实部, 而特征方程 [](#chara) 是 $n$ 次代数方程
$$
\label{algebra_eq}
\lambda ^{n} +a_{1} \lambda ^{n-1} +a_{2} \lambda ^{n-2} +\cdots +a_{n-1} \lambda +a_{n} =0. 
$$

```{prf:theorem} Routh-Hurwitz 准则
由方程 [](#algebra_eq) 的系数作矩阵如下:
$$
\label{matrix}
\left[\begin{array}{cccccc} {a_{1} } & {1} & {0} & {0} & {\cdots } & {0} \\ {a_{3} } & {a_{2} } & {a_{1} } & {1} & {\cdots } & {0} \\ {a_{5} } & {a_{4} } & {a_{3} } & {a_{2} } & {\cdots } & {0} \\ {a_{7} } & {a_{6} } & {a_{5} } & {a_{4} } & {\cdots } & {0} \\ {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } & {\cdots } \\ {a_{2n-1} } & {a_{2n-2} } & {a_{2n-3} } & {\cdots } & {a_{n+1} } & {a_{n} } \end{array}\right], 
$$
其中主对角线上顺次为 $a_{1} ,a_{2} ,\cdots,a_{n}$, 另外, 如果 $i>n$, 则$a_{i} =0$. 那么方程 [](#algebra_eq) 的一切根均有负数实部的充要条件为 [](#matrix) 的各阶主子式大于零, 即:
$$
\Delta _{1} =a_{1} >0,  \Delta _{2} =\left|\begin{array}{cc} {a_{1} } & {1} \\ {a_{3} } & {a_{2} } \end{array}\right|>0, \cdots , \Delta _{n} =a_{n} \Delta _{n-1} >0.
$$
```

对于一个充分光滑的非线性的矢量场 $\dot{x}=f(x)$, 假设 $f(0)=0$, 那么在 $x=0$ 附近, 矢量场具有形式:
$$
\label{nolinear}
\dot{x}=Ax+O(||x||^2),
$$
其中 $A=df(0)$ 为 $f$ 在 $0$ 处的 Jacobi 矩阵. 如果仅在 $0$ 点附近的小邻域内考虑动力学, 线性系统 $\dot{x}=Ax$ 与非线性系统 [](#nolinear) 的零点是否稳定性相同呢? 在此邻域内, 两者的动力学是否是轨道等价的呢?

```{prf:definition}
对于 $C^{1}$ 光滑的矢量场 $\dot{x}=f(x)$, 其平衡点 $x^{*}$ 称为是双曲的, 如果 $df(x^*)$ 的所有特征值都具有非零实部.
```

```{prf:theorem} Hartman-Grobman 定理
对于 $C^{1}$ 光滑的矢量场 $\dot{x}=f(x)$, 如果其平衡点 $x^{*}$ 称为是双曲的, 那么系统在 $x^{*}$ 的一个邻域内拓扑等价于线性系统 $\dot{x}=df(x^*)x$.
```

```{prf:remark}
拓扑等价强于轨道等价, 因此我们可以说, 鞍点, 节点, 和焦点是结构稳定的, 而中心不是结构稳定的. 例如无阻尼线性振子, 阻尼的引入可使其从中心变为焦点.
```