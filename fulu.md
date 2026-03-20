# 附录: 多元微积分


## 隐函数定理

```{prf:theorem}隐函数定理
设 $f : \mathbb{R}^{n+m} \to \mathbb{R}^m$ 为一个连续可微函数, 并设 $\mathbb{R}^{n+m}$ 具有坐标 $(\mathbf{x}, \mathbf{y})$. 固定一点 $(\mathbf{a}, \mathbf{b}) = (a_1, \dots, a_n, b_1, \dots, b_m)$ 使得 $f(\mathbf{a}, \mathbf{b}) = \mathbf{0}$, 其中 $\mathbf{0} \in \mathbb{R}^m$ 是零向量. 

如果雅可比矩阵:

$$J_{f, \mathbf{y}}(\mathbf{a}, \mathbf{b}) = \left[ \frac{\partial f_i}{\partial y_j}(\mathbf{a}, \mathbf{b}) \right]$$

是可逆的, 那么存在一个包含 $\mathbf{a}$ 的开集 $U \subset \mathbb{R}^n$, 使得存在唯一的函数 $g : U \to \mathbb{R}^m$, 满足 $g(\mathbf{a}) = \mathbf{b}$, 且对于所有 $\mathbf{x} \in U$, 都有 $f(\mathbf{x}, g(\mathbf{x})) = \mathbf{0}$. 

此外, $g$ 是连续可微的. 将前一节所示雅可比矩阵的左半部分记为:

$$J_{f, \mathbf{x}}(\mathbf{a}, \mathbf{b}) = \left[ \frac{\partial f_i}{\partial x_j}(\mathbf{a}, \mathbf{b}) \right]$$

则 $g$ 在 $U$ 中的偏导数雅可比矩阵由以下矩阵乘积给出:

$$\left[ \frac{\partial g_i}{\partial x_j}(\mathbf{x}) \right]_{m \times n} = -[J_{f, \mathbf{y}}(\mathbf{x}, g(\mathbf{x}))]_{m \times m}^{-1} [J_{f, \mathbf{x}}(\mathbf{x}, g(\mathbf{x}))]_{m \times n}.$$
```

当 $n=m=1$, 上面的结论就退化为: 如果 $f(x,y)=0$ 是平面上的光滑曲线, 且在某个点 $d_y f(x,y)\neq 0$, 那么这个在这个点附近, 曲线可以表示为一个函数 $y=g(x)$ 的图像, 即 $f(x,g(x))=0$ 在这个点附近局部地成立.

例如球面 $\mathbb{S}^2=\{(x,y,z): ||(x,y,z)||=1\}$ 可以写成是函数 $f(x,y,z)-1=0$, 其中 $f(x,y,z)=||(x,y,z)||$. 当 $z>0$ 时, $d_z f(x,y,z)=2z>0$, 那么由隐函数定理, 球面可以写成是函数图像 $z=g(x,y)$, 实际上, 我们知道 $g(x,y)=\sqrt{1-x^2-y^2}$.


## 多元函数的泰勒展开

函数 $f: \mathbb{R}^n \to \mathbb{R}$ 在 $\mathbf{a} \in \mathbb{R}^n$ 处可微, 当且仅当存在一个线性函数 $L: \mathbb{R}^n \to \mathbb{R}$ 和一个函数 $h: \mathbb{R}^n \to \mathbb{R}$ 满足:

$$f(\mathbf{x}) = f(\mathbf{a}) + L(\mathbf{x} - \mathbf{a}) + h(\mathbf{x})\|\mathbf{x} - \mathbf{a}\|, \quad \lim_{\mathbf{x} \to \mathbf{a}} h(\mathbf{x}) = 0$$

在这种情况下, $L = df(\mathbf{a})$ 是 $f$ 在点 $\mathbf{a}$ 处的微分. 此外, $f$ 在 $\mathbf{a}$ 处的偏导数存在, 且 $f$ 在 $\mathbf{a}$ 处的微分由下式给出:

$$df(\mathbf{a})(\mathbf{v}) = \frac{\partial f}{\partial x_1}(\mathbf{a})v_1 + \dots + \frac{\partial f}{\partial x_n}(\mathbf{a})v_n$$

引入多重指标记法:

$$|\alpha| = \alpha_1 + \dots + \alpha_n, \quad \alpha! = \alpha_1! \dots \alpha_n!, \quad \mathbf{x}^\alpha = x_1^{\alpha_1} \dots x_n^{\alpha_n}$$

对于 $\alpha \in \mathbb{N}^n$ 和 $\mathbf{x} \in \mathbb{R}^n$. 如果 $f: \mathbb{R}^n \to \mathbb{R}$ 的所有 $k$ 阶偏导数在 $\mathbf{a} \in \mathbb{R}^n$ 处连续, 那么根据克莱罗定理, 可以交换 $\mathbf{a}$ 处混合导数的顺序, 因此简写记号:

$$D^\alpha f = \frac{\partial^{|\alpha|} f}{\partial \mathbf{x}^\alpha} = \frac{\partial^{\alpha_1 + \dots + \alpha_n} f}{\partial x_1^{\alpha_1} \dots \partial x_n^{\alpha_n}}$$

在这种情况下的高阶偏导数定义是合理的. 如果 $f$ 的所有 $(k-1)$ 阶偏导数在 $\mathbf{a}$ 的某个邻域内存在且在 $\mathbf{a}$ 处可微, 情况也是如此. 此时我们称 $f$ 在点 $\mathbf{a}$ 处是 $k$ 次可微的.

```{prf:theorem}多元函数泰勒定理
$f: \mathbb{R}^n \to \mathbb{R}$ 在点 $\mathbf{a} \in \mathbb{R}^n$ 处是 $k$ 次连续可微函数. 则存在函数 $h_\alpha: \mathbb{R}^n \to \mathbb{R}$ (其中 $|\alpha| = k$), 使得:

$$
\label{taylor}
f(\mathbf{x}) = \sum_{|\alpha| \leq k} \frac{D^\alpha f(\mathbf{a})}{\alpha!}(\mathbf{x} - \mathbf{a})^\alpha + \sum_{|\alpha|=k} h_\alpha(\mathbf{x})(\mathbf{x} - \mathbf{a})^\alpha,$$

且 $\lim_{\mathbf{x} \to \mathbf{a}} h_\alpha(\mathbf{x}) = 0$. 

如果函数 $f: \mathbb{R}^n \to \mathbb{R}$ 在某个闭球 $B = \{ \mathbf{y} \in \mathbb{R}^n : \|\mathbf{a} - \mathbf{y}\| \leq r \}$ ($r > 0$) 内是 $k + 1$ 次连续可微的, 那么可以根据该邻域内 $f$ 的 $(k+1)$ 阶偏导数推导出余项的精确公式. 即:

$$f(\mathbf{x}) = \sum_{|\alpha| \leq k} \frac{D^\alpha f(\mathbf{a})}{\alpha!}(\mathbf{x} - \mathbf{a})^\alpha + \sum_{|\beta|=k+1} R_\beta(\mathbf{x})(\mathbf{x} - \mathbf{a})^\beta,$$

$$R_\beta(\mathbf{x}) = \frac{|\beta|}{\beta!} \int_0^1 (1 - t)^{|\beta|-1} D^\beta f(\mathbf{a} + t(\mathbf{x} - \mathbf{a})) \, dt.$$

在这种情况下, 由于 $(k+1)$ 阶偏导数在紧集 $B$ 上的连续性, 可以立即得到一致估计:

$$|R_\beta(\mathbf{x})| \leq \frac{1}{\beta!} \max_{|\alpha|=|\beta|} \max_{\mathbf{y} \in B} |D^\alpha f(\mathbf{y})|, \quad \mathbf{x} \in B.$$
```



如果 $\mathbf{f}:\mathbb{R}^{n}\rightarrow\mathbb{R}^m$, $\mathbf{f}$ 的形式如下:
$$
\mathbf{f}(\mathbf(x))=\Big(f_1(\mathbf(x)),f_2(\mathbf(x)),\cdots,f_m(\mathbf(x)) \Big).
$$
那么 $\mathbf{f}$ 展开即只需要对 $m$ 个分量 $f_k$ 作上面的展开就好.

分岔课程中常常把矢量场或者映射的泰勒公式写成多重线性型的形式. 不妨假设 $\mathbf{f}:\mathbb{R}^{n}\rightarrow\mathbb{R}^m$ 是充分光滑的, 那么在 $\mathbf{a}\in \mathbb{R}^{n}$ 处,
$$
\mathbf{f}(\mathbf{x})=\mathbf{a}+A(\mathbf{x}-\mathbf{a})+\frac{1}{2!}B(\mathbf{x}-\mathbf{a},\mathbf{x}-\mathbf{a})+\frac{1}{3!}C(\mathbf{x}-\mathbf{a},\mathbf{x}-\mathbf{a},\mathbf{x}-\mathbf{a})+...
$$
其中对于 $\mathbf{\xi},\mathbf{\eta},\mathbf{\gamma}\in\mathbb{R}^n$, 
$$
A(\mathbf{\xi})=\Big(A_{1}(\mathbf{\xi}),\cdots,A_{m}(\mathbf{\xi})\Big)^T,\\
B(\mathbf{\xi},\mathbf{\eta})=\Big(B_{1}(\mathbf{\xi},\mathbf{\eta}),\cdots,B_{m}(\mathbf{\xi},\mathbf{\eta})\Big)^T,\\
$$
其中 
$$
A_k(\mathbf{\xi})=\sum_{i_1=1}^{1}\frac{\partial f_k}{\partial x_{i_1}}(\mathbf{a})\xi_{i_1},\\
B_k(\mathbf{\xi},\mathbf{\eta})=\sum_{i_1=1}^{n}\sum_{i_2=1}^{n}\frac{\partial ^2f_k}{\partial x_{i_1}\partial x_{i_2}}(\mathbf{a})\xi_{i_1}\eta_{i_2},\\
C_k(\mathbf{\xi},\mathbf{\eta},\mathbf{\gamma})=\sum_{i_1=1}^{n}\sum_{i_2=1}^{n}\sum_{i_3=1}^{n}\frac{\partial ^3f_k}{\partial x_{i_1}\partial x_{i_2}\partial x_{i_3}}(\mathbf{a})\xi_{i_1}\eta_{i_2}\gamma_{i_3}.
$$

如果使用 Einstein 求和约定, 则可以省去求和符号. 例如,
$$
C_k(\mathbf{\xi},\mathbf{\eta},\mathbf{\gamma})=\frac{\partial ^3f_k}{\partial x_{i_1}\partial x_{i_2}\partial x_{i_3}}(\mathbf{a})\xi_{i_1}\eta_{i_2}\gamma_{i_3}.
$$