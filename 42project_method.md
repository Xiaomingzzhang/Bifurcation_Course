# 投影法计算中心流形与约化方程

上一小节我们介绍了中心流形定理, 它允许我们将高维的, 仅有少数特征值退化的方程降阶至低阶方程. 然而, 中心流形的计算是非常繁杂的, 需要进行坐标变换, 以及待定系数法求解中心流形参数表达式的系数. 此前介绍低维系统的分岔时, 我们知道仅仅有几个关键参数控制分岔的方向与类型 (如 Hopf 分岔的 $c_1(0)$), 那么对于一个高维系统, 能否利用中心流形定理, 直接用原系统的一些量表示这些参数, 这样我们可以直接对原系统计算这些参数, 就可判断分岔类型与分岔方向. 这样的方法是存在的, 即我们这节介绍的投影法, 它也是许多分岔软件包进行数值计算分岔类型的基础.

投影法允许我们对一整类系统, 给出中心流形与约化方程的各阶项的系数.

## 线性代数回顾
对于矩阵 $A:\mathbb{R}^{n}\rightarrow\mathbb{R}^n$, 其像空间与零空间分别定义为:
$$
\mathcal{R}(A)&=\{Ax:x\in\mathbb{R}^n\},\\
\mathcal{N}(A)&=\{x:x\in\mathbb{R}^n,Ax=0\}.
$$
注意到, $\mathcal{R}(A),\mathcal{N}(A)$ 均是 $A$ 的不变子空间. 实际上, 也可以认为 $\mathcal{R}(A)$ 是由 $A$ 的非零特征值相应的广义特征矢量所张成, $\mathcal{N}(A)$ 由 $A$ 的零特征值相应的广义特征矢量所张成, 这些子空间作为特征子空间自然是不变的.

```{prf:theorem}
$$
\mathcal{R}(A)=\big(\mathcal{N}(A^{T})\big)^{\perp}.
$$
```
如果 $A$ 是复矩阵, 那么上式中的转置变为共轭转置. 对于子空间 $M$, $M^{\perp}$ 即是 $M$ 的正交补:
$$
M^{\perp}=\{x:\forall y\in M,\langle x,y\rangle=0\}.
$$

例如对于矩阵 $A$, 它有一个一重的零特征值. 我们想找到一个垂直于它的像空间的矢量, 那么根据上述定理, 满足 $A^{T}p=0$ 的非零矢量 $p$ 满足
$$
p\perp \mathcal{R}(A).
$$

在使用投影法计算中心流形的系数时, 需要限制在某个子空间上求逆, 而线性映射本身在这个子空间上是可逆的, 那么如何求逆呢?
```{prf:lemma}
假设矩阵 $A:\mathbb{R}^n\rightarrow\mathbb{R}^n$ 具有不变的子空间 $M$, 且 $A$ 限制在 $M$ 上可逆, 那么对于 $\xi\in M$,
$$
(A|_M)^{-1}\xi=V(V^{T}AV)^{-1}V^{T}\xi,
$$
其中 $V=[v_1,v_2,\cdots,v_k]$, $\{v_1,v_2,\cdots,v_k\}$ 是 $M$ 的一组基.
```

```{prf:proof}
:numbered: false

不妨设 $M$ 由 $\{v_1, \dots, v_k\}$ 张成. 
令 $V = [v_1, \dots, v_k]$ 为 $n \times k$ 的矩阵, 那么:
$$
Av_j = \sum_{i=1}^k B_{ij} v_i,
$$
即 $AV = VB$. $k \times k$ 的矩阵 $B$ 即是 $A|_M$ 在基 $\{v_1, \dots, v_k\}$ 下的矩阵表示. 故 $B$ 是可逆的. 

由 $V^T AV = V^T VB$, 得到
$$
B = (V^T V)^{-1}(V^T AV) \implies B^{-1} = (V^T AV)^{-1} V^T V.
$$

不妨设 $\xi$ 在 $\{v_1, \dots, v_k\}$ 下的表示为 $\hat{\xi}$. 那么
$$
\xi = V\hat{\xi} \implies \hat{\xi} = (V^T V)^{-1} V^T \xi.
$$

故 $B^{-1}\hat{\xi}$ 的表达式为
$$
B^{-1}\hat{\xi} = (V^T AV)^{-1} V^T V \hat{\xi}=(V^T AV)^{-1} V^T \xi.
$$

$B^{-1}\hat{\xi}$ 得到的是 $(A|_M)^{-1}\xi$ 在基 $\{v_1, \dots, v_k\}$ 下的表示. 再换回正常基下的表示, 即 $V B^{-1}\hat{\xi}$. 最终得到:
$$
(A|_M)^{-1}\xi = V(V^T AV)^{-1} V^T \xi.
$$
```
## 单零特征值分岔
首先, 我们作出一些假设, 主要是确保在分岔点处仅有单个 $0$ 特征值.

考虑 $n$ 维含单个参数的自治系统
$$
\dot{z}=f(z,\mu),
$$
其中 $z\in\mathbb{R}^n,\mu\in\mathbb{R}$, $f$ 关于 $z,\mu$ 是充分光滑的矢量场.
```{prf:assumption}
- $f(0,0)=0$;
- $f_z(0,0)$ 具有一个简单的零特征值 (即零特征值的代数重数为 1), 而其他特征值具有负实部.
```

在以上的假设下, 将矢量场关于 $z=0$ 展开为
$$
\dot{z}=\xi(\mu)+A(\mu)z+\frac{1}{2}B_{\mu}(z,z)+\frac{1}{6}C_{\mu}(z,z,z)+\cdots,
$$
其中 $B_\mu,C_\mu$ 分别为 $f$ 在 $z=0$ 处的二次和三次线性型, 下标表示其对参数 $\mu$ 的依赖性. 将 $A(\mu),\xi(\mu)$ 进一步写为
$$
A(\mu)=&A_0+A_1\mu+O(\mu^2),\quad A_0=A(0),A_1=A'(0),\\
\xi(\mu)=&\xi_1\mu+O(\mu^2),\quad \xi_1=\xi'(0).
$$

下面我们将对矢量场的状态变量分别投影到零特征值相应的特征矢量方向, 即 $\mathcal{N}(A_0)$, 和其他具有负实部的特征值对应的特征空间, 即 $\mathcal{R}(A_0)$. 

具体地, 取矢量 $p_0,q_0\in\mathbb{R}^n$ 满足:
$$
A_0q_0=0,A_0^{T}p_0=0,\langle q_0,p_0\rangle=1.
$$
那么由于 $p_0\in\mathcal{N}(A^{T})$, $p_0\perp \mathcal{R}(A_0)$.

对于任意矢量 $z\in\mathbb{R}^n$, 作如下分解:
$$
z=xq_0+y,x\in\mathbb{R},y\in\mathcal{R}(A_0)=\text{span}\{p_0\}^{\perp}.
$$
如图所示:
```{figure} ./asserts/figs/投影法示意图.png
:alt: 图片无法加载
:width: 80%
:align: center
```
那么我们得到:
$$
x=\langle p_0,z\rangle,y=z-\langle p_0,z\rangle q_0.
$$

实际上, 如果我们定义线性映射:
$$
P_n z=\langle p_0,z\rangle q,P_s z=z-\langle p_0,z\rangle q_0.
$$
可以得到:
$$
P_n^2=P_n,P_s^2=P_z.
$$
故 $P_n,P_z$ 都为投影算子, 分别将 $z$ 投影到 $q_0$ 张成的子空间与 $\mathcal{R}(A_0)$ 上. 那么可以将 $x,y$ 的表达式重写为
$$
P_nz=xq_0,P_sz=y.
$$
```{prf:remark}
一个方阵 $P$ 满足 $P^2=P$, 则称其为一个投影算子, 如果它还满足 $P^T=P$, 那么这个投影是正交的. 实际上, 此时 $\mathcal{N}(P)^{\perp}=\mathcal{R}(P)$, 也就是零空间与像空间互为正交补.
```
关于 $x,y$ 求导得到:
$$
\dot{x}&=\langle p_0,f(xq_0+y,\mu)\rangle,\\
\dot{y}&=f(xq_0+y,\mu)-\langle p_0,f(xq_0+y,\mu)\rangle q_0.
$$
注意到, 上面的方程是 $n+1$ 维的方程, 但是由于 $y$ 自然地约束在 $n-1$ 维的子空间上, 本质上并没有增加额外的维度.

代入 $f$ 的多重线性型展开式, 并略去与范式的无关项, 我们得到:
$$
\label{pm_danling}
\dot{x}&=\alpha\mu+a\mu x+\frac{1}{2}cx^2+\langle p_0,B(q_0,y)\rangle x+\frac{1}{6}\delta x^3+\cdots,\\
\dot{y}&=A_0 y+\frac{1}{2}H_{20}x^2+\cdots,
$$
其中
$$
\alpha&=\langle p_0,\xi_1\rangle,\\
a&=\langle p_0,A_1 q_0\rangle,\\
c&=\langle p_0,B( q_0,q_0)\rangle,\\
\delta&=\langle p_0,C(q_0, q_0,q_0)\rangle,\\
H_{20}&=B( q_0,q_0)-\langle p_0,B( q_0,q_0)\rangle q_0,
$$
其中多重线性型 $B=B_0,C=C_0$ 为其在 $\mu=0$ 处的取值.
```{prf:remark}
保留以上推导各项的原则是保证约化方程里的项 $\mu,\mu x,x^2,x^3$ 的系数是已知的. 由于后面要代入中心流形 $y=h(x,\mu)=O(|(x,\mu)|^2)$, 利用这一点也可舍去一些不需要的项.
```
由 [](#pm_danling) 的第一式, 中心流形 $y=h(x,\mu)$ 中, 我们仅需要计算 $x^2$ 项的系数即可. 假设中心流形的表达式为
$$
y=h(x,\mu)=\frac{1}{2}h_{20}x^2+\cdots.
$$
那么由不变方程:
$$
\dot{y}=h_x(x,\mu)\dot{x},
$$
我们得到
$$
A_0 (\frac{1}{2}h_{20}x^2+\cdots)+\frac{1}{2}H_{20}x^2+\cdots=(h_{20}x+\cdots)\big(\alpha\mu+a\mu x+\frac{1}{2}cx^2+\langle p_0,B(q_0,h(x,\mu))\rangle x+\frac{1}{2}\delta x^3+\cdots\big),
$$
那么有
$$
A_0h_{20}=-H_{20}.
$$
由 $H_{20}$ 的表达式, $H_{20}\in\mathcal{R}(A_0)$. 并且由于 $y=h(x,\mu)\in\mathcal{R}(A_0)$, 这表明 $h(x,\mu)$ 中的系数向量也必须在子空间 $\mathcal{R}(A_0)$, 故只需限制在 $\mathcal{R}(A_0)$ 上求逆即可得到 $h_{20}$.

因此中心流形为
$$
h(x,\mu)=-\frac{1}{2}A_0^{-1}H_{20}x^2+\cdots.
$$
代入到 $x$ 的微分方程中, 我们得到
$$
\dot{x}=\alpha\mu+a\mu x+\frac{1}{2}cx^2+\frac{1}{6}d x^3+\cdots,
$$
其中
$$
\alpha&=\langle p_0,\xi_1\rangle,\\
a&=\langle p_0,A_1 q_0\rangle,\\
c&=\langle p_0,B( q_0,q_0)\rangle,\\
H_{20}&=B( q_0,q_0)-\langle p_0,B( q_0,q_0)\rangle q_0,\\
d&=\langle p_0,C(q_0, q_0,q_0)\rangle-3\langle p_0,B(q_0, A_0^{-1}H_{20})\rangle.
$$

根据此前介绍的结果: [](./21anjie.md), [](./22yincha), 即可判断上述系数判定系统的分岔. 具体实例可参考 [](./61symbolic.md).

## Hopf 分岔

考虑 $n$ 维含单个参数的自治系统
$$
\dot{x}=f(x,\mu),
$$
其中 $x\in\mathbb{R}^n,\mu\in\mathbb{R}$, $f$ 关于 $x,\mu$ 是充分光滑的矢量场.
```{prf:assumption}
- $f(0,\mu)=0$;
- $f_z(0,0)$ 具有一对简单的共轭复特征值:
$$
\lambda(\mu)=\alpha(\mu)+ i\omega(\mu),
$$
且
$$
\alpha(0)=0,\alpha'(0)\neq0,\omega_0=\omega(0)>0,
$$
$f_z(0,0)$ 的其余特征值均具有负实部.
```
