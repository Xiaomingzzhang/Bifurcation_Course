# 周期解与 Floquet 理论

对于自治系统或者非自治的时间周期 $T$ 系统, 一个周期为 $T$ 的周期解是流映射 $\varphi^T$ 的不动点. 判断这个不动点的稳定性需要得到 $\varphi^T$ 在该不动点处的导数, 而这个导数是一个线性周期系数的矩阵微分方程的解, 即变分方程的解. 特征乘子即上述导数的特征值. 对于自治系统, 特征乘子中总存在 1, 而其相应的特征向量为矢量场在该不动点处的取值, 这是由于周期轨道上的所有点都是流映射 $\varphi^T$ 的不动点; 如果在周期轨道附近取 Poincaré 截面, 那么相应的 Poincaré 映射的导数的特征值将会保留其他除了 1 以外的特征值. 上述想法可以自然推广到分段光滑系统和碰撞系统.

## 时间推进映射的导数与变分方程
考虑 $n$ 阶非自治系统:
$$
\label{floquet_basic_eq}
\dot{x}=f(x,t),
$$
其中 $x\in\mathbb R^n$, $f:\mathbb{R}^n\times\mathbb{R}\rightarrow\mathbb{R}^n$ 是 $C^{1}$ 类的矢量场. 我们假设系统的任意解可以延拓到整个实轴. 给定初值 $x(t_0)=x_0$, 记 [](#floquet_basic_eq) 的解为 $\phi(t,x_0,t_0)$. 

**给定时间 $t_0,t_1$**, 定义时间推进映射:
$$
\psi:x_0\mapsto \phi(t_1,x_0,t_0).
$$
记 $\psi$ 的导数为 $\Phi(t_1,t_0)$, 那么
$$
\Phi(t_1,t_0)(x_0)=d\psi(x_0)=\phi_{x_0}(t_1,x_0,t_0).
$$
```{prf:example}
考虑微分方程:
$$
\dot{x}&=-x,\\
\dot{y}&=xy
$$
方程以初值 $x(0)=x_0,y(0)=y_0$ 的解为:
$$
\phi(t,(x_0,y_0),0)=\big(e^{-t}x_0,e^{x_0-e^{-t}x_0}y_0\big).
$$
因此
$$
\Phi(t,0)((x_0,y_0))=\begin{bmatrix}
    e^{-t}&0\\
    (e^t-1)e^{x_0-t-e^{-t}x_0}y_0&e^{x_0-e^{-t}x_0}
\end{bmatrix}.
$$
```
对于自治系统, 假设周期解上一点为 $x^*$, 那么 $\Phi(T,0)(x^*)$ 即时间推进映射的导数, 其中 $T$ 为周期解的周期.

```{prf:theorem}
给定时间 $t_1,t_0$, 系统 [](#floquet_basic_eq) 时间推进映射在 $x_0$ 处的导数 $\Phi(t_1,t_0)(x_0)$ 是以下矩阵微分的初值问题的解在时间 $t_1$ 处的取值:
$$
\label{floquet_var}
\dot{X}=f_x(\phi(t,x_0,t_0),t)X,\quad X(t_0)=I,
$$
其中 $X$ 为 $n\times n$ 的矩阵, $I$ 为单位矩阵.
```
```{prf:proof}
:numbered: false
由 $\phi(t,x_0,t_0)$ 满足 [](#floquet_basic_eq), 代入方程中并两边对于 $x_0$ 求导得到:
$$
\dot{\phi}_{x_0}(t,x_0,t_0)=\frac{d}{dt}\phi_{x_0}(t,x_0,t_0)=f_x(\phi(t,x_0,t_0),t)\phi_{x_0}(t,x_0,t_0).
$$
由 $\phi(t_0,x_0,t_0)=x_0$, 我们有 $\phi_{x_0}(t_0,x_0,t_0)=I$. 故 $\phi_{x_0}(t,x_0,t_0)=\Phi(t,t_0)(x_0)$ 是 [](#floquet_var) 的解.
```
矩阵微分方程 [](#floquet_var) 称作是 [](#floquet_basic_eq) 的解 $\phi(t,x_0,t_0)$ 的变分方程.
```{prf:example}
对于线性时变系统 $\dot{x}=A(t)x$, 其中矩阵 $A(t)$ 的每个分量都关于时间 $t$ 连续, 系统任意一个解的变分方程为:
$$
\dot{X}=A(t)X,
$$
其中 $X$ 为 $n\times n$ 的矩阵.
```

```{prf:proposition}
- $\Phi(t,t_0)=\Phi(t,t_1)\Phi(t_1,t_0)$;
- 若矢量场 $f(x,t)=f(x)$ 不依赖于时间, 则 $\Phi(t,t_0)=\Phi(t-t_0,0)$, 且有
$$
f(\phi(t,x_0,0))=\Phi(t,0)(x_0)f(x_0).
$$
特别地, 如果 $\phi(T,x_0,0)=x_0$, 那么 $\Phi(T,0)(x_0)f(x_0)=f(x_0)$, 即 $f(x_0)$ 是 $\Phi(T,0)(x_0)$ 的 1 特征值相应的特征向量 (前提是 $f(x_0)\neq 0$).
```
```{prf:proof}
:numbered: false
利用非自治系统解的性质:
$$
\phi(t,x_0,t_0)=\phi(t,\phi(t_1,x_0,t_0),t_1),
$$
两边关于 $x_0$ 求导即得到了 $\Phi(t,t_0)(x_0)=\Phi(t,t_1)(\phi(t_1,x_0,t_0))\Phi(t_1,t_0)(x_0)$.

如果 $f$ 不依赖于时间, 那么 $\phi(t,x_0,t_0)=\phi(t-t_0,x_0,0)$, 两边关于 $x_0$ 求导即得到 $\Phi(t,t_0)=\Phi(t-t_0,0)$.

由于 $\phi(t,\phi(s,x_0,0),0)=\phi(s,\phi(s,x_0,0),0)$, 两边关于 $t$ 在 $t=0$ 处求导得到:
$$
\dot{\phi}(t,\phi(s,x_0,0),0)|_{t=0}&=f(\phi(s,x_0,0))\\
&=\phi_x(s,\phi(s,x_0,0),0)\dot{\phi}(t,x_0,0)|_{t=0}\\
&=\phi_x(s,\phi(s,x_0,0),0)f(x_0),
$$
即
$$
\Phi(s,0)(x_0)f(x_0)=f(\phi(s,x_0,0)).
$$
```
## Floquet 定理与特征乘子
由前一小节, 周期解处的变分方程 $\dot{X}=A(t)X$ 的系数矩阵也是周期的, 其中对于自治系统, $A(t)=f_x(\phi(t,x_0,0))$, 对于非自治周期系统, $A(t)=f_x(\phi(t,x_0,t_0),t)$. 本小节我们考虑周期矩阵方程:
$$
\label{floquet_matrix_eq}
\dot{X}=A(t)X,
$$
其中 $X,A$ 是 $n\times n$ 的矩阵, $A$ 的每个分量都是关于时间的周期 $T$ 函数.
```{prf:definition}
若方程 [](#floquet_matrix_eq) 初值的行列式不等于 0, 称其解为基解矩阵; 以 $X(t_0)=I$ 为初值的解 $X(t)$ 称为主解矩阵.
```
注意到, 如果 $X(t)$ 是方程 [](#floquet_matrix_eq) 的解, 那么它的每一列都是方程 $\dot{x}=A(t)x$ 的解.

若存在某个 $t_0$ 使得 $X(t_0)$ 的行列式为 $0$, 那么 $X(t)$ 的行列式永远为 $0$. 实际上, 由假设存在 $c\in\mathbb{R}^n,c\neq 0$ 使得 $X(t_0)c=0$. 由于 $X(t)c$ 也是 $\dot{x}=A(t)x$ 的解, $x(t)\equiv 0$ 与 $X(t)c$ 在时刻 $t_0$ 相交, 由解的存在唯一性, $X(t)c\equiv 0$. 因此, **若 $X(t)$ 在某个时刻行列式不为 0, 那么其行列式永远也不为 0.**

```{prf:lemma}
假设 $n\times n$ 的矩阵 $C$ 非奇异, 那么存在矩阵 $B$ 使得 $e^{B}=C$.
```

```{prf:theorem} Floquet 定理
:label: thm_floquet
[](#floquet_matrix_eq) 的每个基解矩阵的形式如下:
$$
X(t)=P(t)e^{Bt},
$$
其中 $P$ 也是周期 $T$ 的矩阵, $B$ 是常数矩阵.
```
```{prf:proof}
:numbered: false
若 $X(t)$ 是基解矩阵, 那么由于 $A(t)$ 的周期性, $X(t+T)$ 也是基解矩阵. 那么存在非奇异矩阵 $C$ 使得
$$
X(t+T)=X(t)C.
$$
实际上, 取 $C=X(0)^{-1}X(T)$, 那么 $X(t)C$ 与 $X(t+T)$ 在 $t=0$ 处具有相同的初值, 由解的存在唯一性, 上式对于任意的时间 $t$ 恒成立.

由之前的引理, 存在 $B$ 使得 $e^{BT}=C$. 令 $P(t)=X(t)e^{-Bt}$. 那么 $P(t)$ 满足定理的条件, 只需说明 $P(t)$ 是 $T$ 周期的即可.
$$
P(t+T)=X(t+T)e^{-B(t+T)}=X(t)e^{BT}e^{-B(t+T)}=P(t).
$$
```

```{prf:corollary}
存在非奇异的依赖于时间的变换, 将线性方程 $\dot{x}=A(t)x$ 变为常系数方程.
```
```{prf:proof}
:numbered: false
取 $P(t),B$ 如 [](#thm_floquet) 所述. 取坐标变换 $x=P(t)y$. 由于
$$
0=\frac{d}{dt}(P(t)^{-1}P(t))=(\frac{d}{dt}P(t)^{-1})P(t)+P(t)^{-1}\dot{P}(t),
$$
我们有:
$$
\frac{d}{dt}P(t)^{-1}=-P(t)^{-1}\dot{P}(t)P(t)^{-1}.
$$
那么
$$
\dot{y}=P^{-1}(AP-\dot{P})y.
$$
由于 $P(t)=X(t)e^{-Bt}$, 知 $\dot{P}=AP-PB$, 代入上式即得到了结论.
```
```{prf:definition}
如果 $X(t)$ 是 [](#floquet_matrix_eq) 的基解矩阵, 满足 $X(t+T)=X(t)C$ 的矩阵 $C$ 称为单值矩阵. 特别地, 如果 $X(0)=I$, 则 $C=X(T)$. $C$ 的特征值称为是 [](#floquet_matrix_eq) 的特征乘子.
```
```{prf:remark}
特征乘子不依赖于基解矩阵的选择. 实际上, 若 $Y(t)$ 也是基解矩阵, 那么存在非奇异矩阵 $D$ 使得 $Y(t)=X(t)D$ (利用解的存在唯一性). 于是
$$
Y(t+T)=X(t+T)D=X(t)CD=Y(t)D^{-1}CD.
$$
这说明 $Y(t)$ 与 $X(t)$ 的单值矩阵是相似的.
```

```{prf:theorem}
对于 $C^{1}$ 类的矢量场 $\dot{x}=f(x)$, 如果系统存在周期轨道 $x^*(t)$, 并且线性周期系统:
$$
\dot{X}=A(t)X
$$
的特征乘子的模都小于 1 (排除一个为 1 的特征乘子), 那么周期轨道是渐进稳定的, 其中 $A(t)=f_x(x^*(t))$.
```
## 周期轨道附近的 Poincaré 映射
对于自治系统 $\dot{x}=f(x)$, 可在周期轨道附近取一个超平面 $\Sigma$, 使得周期轨道横截地穿过该平面. 在横截点附近的一个邻域 $U$ 内 (限制在该超平面上), 建立 Poincaré 映射:
$$
P:U\mapsto\Sigma,\quad P(x)=\phi(\tau(x),x,0),
$$
其中 $\tau(x)$ 是再次回到 $U$ 的邻域的时间. 实际上, 如果超平面由 $\{x:H(x)=0\}$ 给出, 其中 $H$ 是 $\mathbb{R}^n$ 上的可微函数. 在周期点 $x^*\in U$处,
$$
H(\phi(T,x^*,0))=0,\quad \frac{d}{dt}H(\phi(t,x,0))|_{t=T,x=x^*}=\nabla H(x^*)\cdot f(x^*)\neq 0,
$$
其中 $T$ 是周期解的周期, $\nabla H(x^*)\cdot f(x^*)\neq 0$ 是横截性条件.
由隐函数定理, 存在可微函数 $\tau(x)$ 使得:
$$
\tau(x^*)=T,\quad H(\phi(\tau(x),x,0))=0.
$$

```{prf:theorem}
$\phi_x(T,x^*,0)$ 除 1 以外的特征值由 $dP(x^*)$ 给出. 并且两个矩阵的特征多项式有如下的联系:
$$
\text{det}(\lambda I-\phi_x(T,x^*,0))=(\lambda-1)\text{det}(\lambda I-dP(x^*)).
$$
```
```{prf:proof}
:numbered: false
由映射 $P$ 的定义, 我们有:
$$
\label{floquet_eq_poincare}
dP(x^*)&=\dot{\phi}(\tau(x),x,0)d\tau(x)|_{x=x^*}+\phi_x(\tau(x),x,0)|_{x=x^*}\\
&=f(x^*)d\tau(x^*)+\phi_x(T,x^*,0).
$$
任意矢量 $\eta\in\mathbb{R}^n$ 我们可以将其分解为:
$$
\eta=\xi+af(x^*),
$$
其中 $\xi\in\Sigma,a\in\mathbb{R}$. 对于 $\xi\in\Sigma$, 我们有:
$$
\phi_x(T,x^*,0)\xi=dP(x^*)\xi-(d\tau(x^*)\xi)f(x^*),
$$
其中我们用到了 [](#floquet_eq_poincare) 以及 $d\tau(x^*)\xi$ 为实数. 又由于:
$$
\phi_x(T,x^*,0)f(x^*)=f(x^*).
$$
那么在上述空间分解 $(a,\xi)$ 下, 线性映射 $\phi_x(T,x^*,0)$ 的矩阵具有形式:
$$
\phi_x(T,x^*,0)=
\left(
\begin{array}{c|c}
1 & -d\tau(x^*) \\
\hline
0 & \\
\vdots & dP(x^*) \\
0 &
\end{array}
\right).
$$
```
```{prf:remark}
对于非自治周期系统 $\dot{x}=f(x,t)$, 将其扩充为 $n+1$ 阶的系统, 那么取特定的相位为 Poincaré 截面, 此时 Poincaré 映射去除掉的 1 乘子对应的特征向量就是额外增加的维度的方向.
```

<!-- ## 高阶变分方程

研究周期解的分岔问题, 知道时间推进映射在周期解处的一阶导数往往是不够的, 而计算其高阶导数则需要计算高阶变分方程的解. -->

## Floquet 理论在分段光滑系统与碰撞系统的推广
如果系统存在切换行为或者类似碰撞的现象, 并且系统的轨线横截地穿过切换面或者碰撞面, Floquet 理论仍然是适用的, 不过此时时间推进映射的导数要乘以一个与切换或者碰撞相关的矩阵, 文献中常常称作是跳跃矩阵.
### 分段光滑系统
考虑如下分段光滑系统:
$$
\dot{x}&=f_1(x,t), \quad H(x)<0\\
\dot{x}&=f_2(x,t), \quad H(x)>0
$$
其中 $f_1,f_2,H$ 都至少是 $C^{1}$ 类的, $H<0,H>0$ 时的解分别记为 $\phi_1(t,x_0,t_0),\phi_2(t,x_1,t_1)$. 

假设分段系统的一个解满足:
- $H(\phi_1(t,x^*,t_0))<0$, $t\in[t_0,\tilde{t})$;
- $H(\phi_1(\tilde{t},x^*,t_0))=0$;
- $H(\phi_2(t,\phi_1(\tilde{t},x^*,t_0),\tilde{t}))>0$, $t\in(\tilde{t},t_1]$.

并且解穿过 $H(x)=0$ 时是横截的, 即:
$$
\nabla H(\phi_1(\tilde{t},x^*,t_0))\cdot f(\phi_1(\tilde{t},x^*,t_0))\neq 0.
$$
由先前的论述, 穿越时间 $\tau(x)$ 在 $x^*$ 的一个邻域 $U$ 内是有定义的且是光滑的.

定义时间推进映射:
$$
\psi:U\rightarrow\mathbb{R}^n, \quad x\mapsto \phi_2(t_1,\phi_1(\tau(x),x,t_0),\tau(x)).
$$
```{prf:theorem}
:label: piecewise_floquet
$$
d\psi(x)=\Phi_2(t_1,\tau(x))S\Phi_1(\tau(x),t_0),
$$
其中 $\Phi_1,\Phi_2$ 分别是 $f_1,f_2$ 的基解矩阵, $S$ 称为跳跃矩阵, 由如下公式给出:
$$
S=I+\frac{(f_2-f_1)\nabla H^{T}}{\langle \nabla H,f_1\rangle},
$$
其中 $f_1,f_2,\nabla H$ 均在解接触到 $H(x)=0$ 处取值.
```
```{prf:remark}
- 上面的定理中, $\Phi_1(\tau(x),t_0)$ 在 $x$ 处取值, $\Phi_2(t_1,\tau(x))$ 在 $\phi_1(\tau(x),x,t_0)$ 处取值;
- 可以看到, 对于分段系统, 基解矩阵的表达式只不过是在两个矢量场切换之间加了一个跳跃矩阵 $S$. 特别地, 当 $f_2=f_1$ 时, $S$ 为单位矩阵;
- 如果 $\nabla H\cdot f_2\neq 0$, 那么 $S$ 是非退化的, 且:
$$
S^{-1}=I+\frac{(f_1-f_2)\nabla H^T}{\langle \nabla H,f_2\rangle}.
$$
```
```{prf:example}
考虑分段光滑系统
$$
\dot{x}&=y\\
\dot{y}&=-kx+f(t),
$$
其中,
- 当 $x<c$ 时, $k=k_1$;
- 当 $x>c$ 时, $k=k_2$.

自然地, $H(x,y)=x-c$. 那么当解从 $x<c$ 穿过 $x=c$ 到达 $x>c$ 时, 跳跃矩阵为:

$$
S=
\begin{bmatrix}
    1&0\\
    0&1
\end{bmatrix}
+\frac{1}{y^*}
\begin{bmatrix}
    0&0\\
    c(-k_2+k_1)&0
\end{bmatrix},
$$
其中 $y^*$ 为解穿过 $x=c$ 时的 $y$ 的取值.
```
在下面的引理的帮助下, [](#piecewise_floquet) 的证明是比较简单的.
```{prf:lemma}
在上面的记号下, 我们有:
$$
d\tau(x)=-\frac{1}{\nabla H\cdot f_1}\nabla H^{T}\Phi_1(\tau(x),t_0),\\
\frac{\partial }{\partial t_0}\phi_1(t,x,t_0)=-\big(\frac{\partial }{\partial x}\phi_1(t,x,t_0)\big)f_1(x,t_0).
$$
```
```{prf:proof}
:numbered: false
由于
$$
H(\phi_1(\tau(x),x,t_0))=0,
$$
关于 $x$ 求导即得到第一个等式.

关于等式
$$
\phi_1(t,x,t_0)=\phi_1(t,\phi_1(t_0+s,x,t_0),t_0+s)
$$
在 $s=0$ 处求导得到:
$$
0=&\phi_{1x}(t,\phi_1(t_0+s,x,t_0),t_0+s)\dot{\phi}_1(t_0+s,x,t_0)|_{s=0}+\phi_{1t_0}(t,\phi_1(t_0+s,x,t_0),t_0+s)|_{s=0}\\
&=\phi_{1x}(t,x,t_0)f(x,t_0)+\phi_{1t_0}(t,x,t_0).
$$
```
```{prf:proof} [](#piecewise_floquet) 的证明
:numbered: false
$$
\frac{d}{dx}\psi(x)&=\frac{d}{dx}\phi_2(t_1,\phi_1(\tau(x),x,t_0),\tau(x))\\
&=\phi_{2x}(t_1,\phi_1(\tau(x),x,t_0),\tau(x))\big(\dot{\phi}_1(\tau(x),x,t_0)d\tau(x)+\phi_{1x}(\tau(x),x,t_0)\big)\\
&+\phi_{2t_0}(t_1,\phi_1(\tau(x),x,t_0),\tau(x))d\tau(x)\\
&=\Phi_2(t_1,\tau(x))\big(f_1d\tau(x)+\Phi_1(\tau(x),t_0)\big)-\Phi_2(t_1,\tau(x))f_2d\tau(x)\\
&=\Phi_2(t_1,\tau(x))\Phi_1(\tau(x),t_0)+\Phi_2(t_1,\tau(x))(f_1-f_2)d\tau(x)\\
&=\Phi_2(t_1,\tau(x))S\Phi_1(\tau(x),t_0).
$$
```
```{prf:remark}
假设解在时间 $t_i,t_{i+1}$ 位于第 $k_i$ 个矢量场 ($i=1,\cdots,m$), 并且每次穿越都是横截的, 那么我们有:
$$
d\psi(x)=\Phi_{k_i}(t_{i+1},t_{i})S_{i}\Phi_{k_{i-1}}(t_{i},t_{i-1})S_{i-1}\Phi_{k_{i-2}}(t_{i-1},t_{i-2})\cdots \Phi_{k_1}(t_{2},t_{1}).
$$
```
[](#piecewise_floquet) 对于可能出现滑移的情况也是适用的. 如果在分界面两侧的矢量场具有 "挤向" 分界面的趋势, 我们容许轨线可以在分界面上进行移动. 现实中, 这对应于带有干摩擦的系统. 
```{prf:proposition}
如果 $f_2$ 定义在超曲面 $\{x:H(x)=0\}$ 上, 那么跳跃矩阵 $S$ 的零空间由 $f_2-f_1$ 张成.
```
```{prf:proof}
:numbered: false
由于 $\nabla H\cdot f_2=0$, 我们有
$$
S(f_2-f_1)=f_2-f_1+(f_2-f_1)\frac{\nabla H^{T}(f_2-f_1)}{\langle \nabla H,f_1\rangle}=0.
$$
而一个形如 $I+ab^T$ 的矩阵的特征值为 $\{1+b^Ta,1,\cdots,1\}$, 其中 $a,b$ 都为列向量. 令:
$$
a=\frac{f_2-f_1}{\langle \nabla H,f_1\rangle},\quad b=\nabla H.
$$
那么 $b^Ta=-1$. 故 $S$ 仅有一个单重的 $0$ 特征值.
```
```{prf:remark}
对于矩阵 $I+ab^T$ ($a,b\neq 0$), 与 $b$ 正交的子空间的基 $\{u_1,\cdots,u_{n-1}\}$ 构成了其 $n-1$ 个 $1$ 相应的特征向量. 又由于:
$$
(I+ab^T)a=a+ab^Ta=a(1+b^Ta)=(1+b^Ta)a,
$$
故 $a$ 是特征值 $1+b^Ta$ 相应的特征向量.
```
上面的命题说明对于干摩擦系统, 如果出现周期轨道出现滑移段, 那么周期解总存在一个 0 特征乘子.
### 碰撞系统

考虑系统 $\dot{x}=f(x,t)$. 假设系统的相空间中存在一个超曲面 $\{x:H(x)=0\}$, 以及一个保持该曲面不变的映射 $G$, 即 $y\in\{x:H(x)=0\}$, 则有 $H(G(y))=0$.

系统的轨线到达曲面 $\{x:H(x)=0\}$ 后, 其状态发生如下变化:
$$
(t,x)\mapsto (t,G(x)).
$$

同样假设轨线从曲面一侧出发, 横截地碰到曲面, 经过状态变化后, 由经历一段时间在相空间中移动. 设开始与结束的时间分别为 $t_0,t_1$, 初始点为 $x^*$, 那么由横截性条件, 在 $x^*$ 附近存在一个邻域 $U$ 使得到达时间 $\tau(x)$ 在此邻域内定义良好.

定义时间推进映射:
$$
\psi:U\rightarrow \mathbb{R}^n,\quad \psi(x)=\phi(t_1,G(\phi(\tau(x),x,t_0)),\tau(x)).
$$
```{prf:theorem}
$$
d\psi=\Phi(t_1,\tau(x))S\Phi(\tau(x),t_0),
$$
其中
$$
S=dG(x_1)+\frac{\big(f(G(x_1),\tau(x))-dG(x_1)f(x_1,\tau(x))\big)\nabla H(x_1)^T}{\langle \nabla H(x_1),f(x_1,\tau(x))\rangle}
$$
其中 $x_1$ 为轨线刚刚到达 $\{x:H(x)=0\}$ 的状态.
```
```{prf:proof}
:numbered: false
$$
\frac{d}{dx}\psi(x)&=\frac{d}{dx}\phi(t_1,G(\phi(\tau(x),x,t_0)),\tau(x))\\
&=\phi_x(t_1,G(x_1),\tau(x))dG(x_1)\big(\dot{\phi}(\tau(x),x,t_0)d\tau(x)+\phi_x(\tau(x),x,t_0)\big)
+\phi_{t_0}(t_1,G(x_1),\tau(x))d\tau(x)\\
&=\Phi(t_1,\tau(x))dG(x_1)\big(f(x_1,\tau(x))d\tau(x)+\Phi(\tau(x),t_0)\big)
-\Phi(t_1,\tau(x))f(G(x_1),\tau(x))\\
&=\Phi(t_1,\tau(x))dG(x_1)\Phi(\tau(x),t_0)+\Phi(t_1,\tau(x))\big(dG(x_1)f(x_1,\tau(x))-f(G(x_1),\tau(x))\big)d\tau(x)\\
&=\Phi(t_1,\tau(x))S\Phi(\tau(x),t_0).
$$
```
```{prf:example}
考虑单自由度碰撞振子:
$$
\dot{x}&=y
\dot{y}&=-kx-cy+f(t)
$$
$H(x,y)=x, G(x,y)=(x,-ry)^{T}$, 其中 $0<r\leq 1$ 为恢复系数. 设碰撞状态为 $(y^*,t^*)$, 则

$$
dG=\begin{bmatrix}
    1&0\\
    0&-r
\end{bmatrix},
\nabla H=\begin{bmatrix}
    1\\
    0
\end{bmatrix},
f(0,y^*,t^*)=
\begin{bmatrix}
    y^*\\
    -cy^*+f(t^*)
\end{bmatrix},
f(0,-ry^*,t^*)=
\begin{bmatrix}
    y^*\\
    cry^*+f(t^*)
\end{bmatrix}.
$$
此时有

$$
S=\begin{bmatrix}
    -r&0\\
    \frac{f(t^*)+rf(t^*)}{y^*}&-r
\end{bmatrix}.
$$
```