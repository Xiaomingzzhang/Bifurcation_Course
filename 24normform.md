

# 范式方法

## 微分方程: $(ad_A)^*=ad_{A^*}$ 的证明
为了证明线性算子 $ad_A$ 在所定义的内积下的共轭为 $ad_{A^*}$, 我们需要证明对于任意向量值齐次多项式 $P(x), Q(x) \in H_n^k$, 以下等式成立:

$$\langle ad_A(P), Q\rangle = \langle P, ad_{A^*}(Q)\rangle$$

这里, $P(x) = (p_1(x), \dots, p_n(x))^T$ 且 $Q(x) = (q_1(x), \dots, q_n(x))^T$. 向量场的内积按分量定义为 $\langle P, Q\rangle = \sum_{i=1}^n \langle p_i, q_i\rangle$, 其中标量内积为 $\langle f, g\rangle = f(\partial)g(x)|_{x=0}$.

我们首先给出一个重要性质: $x_m$ 和偏导数 $\frac{\partial}{\partial x_m}$ 在该内积下是共轭算子. 对于任意标量多项式 $f$ 和 $g$:

* **性质 A:** $\langle x_m f, g\rangle = \langle f, \frac{\partial g}{\partial x_m}\rangle$
    *证明:* $\langle x_m f, g\rangle = \frac{\partial}{\partial x_m} f(\partial) g(x)|_{x=0} = f(\partial) (\frac{\partial g}{\partial x_m})|_{x=0} = \langle f, \frac{\partial g}{\partial x_m}\rangle$.
* **性质 B:** $\langle \frac{\partial f}{\partial x_m}, g\rangle = \langle f, x_m g\rangle$
    *证明:* 由于内积的对称性, $\langle \frac{\partial f}{\partial x_m}, g\rangle = \langle g, \frac{\partial f}{\partial x_m}\rangle = \langle x_m g, f\rangle = \langle f, x_m g\rangle$.

根据定义, $ad_A(P(x)) = DP(x)Ax - AP(x)$. 
我们可以将其分解为两个线性算子: $L_1(P) = DP(x)Ax$ 和 $L_2(P) = AP(x)$. 我们将分别求出两者的共轭.

### 求解 $L_2(P) = AP(x)$ 的共轭
设矩阵 $A$ 的元素为 $a_{ij}$. $AP(x)$ 的第 $i$ 个分量为 $\sum_{j=1}^n a_{ij} p_j(x)$.

$$\langle L_2(P), Q\rangle = \langle AP, Q\rangle = \sum_{i=1}^n \langle \sum_{j=1}^n a_{ij} p_j, q_i\rangle = \sum_{i=1}^n \sum_{j=1}^n a_{ij} \langle p_j, q_i\rangle$$

通过交换求和顺序, 我们得到:

$$\langle L_2(P), Q\rangle = \sum_{j=1}^n \langle p_j, \sum_{i=1}^n a_{ij} q_i\rangle$$

假设为实内积空间, 转置矩阵 $A^*$ 的元素为 $(A^*)_{ji} = a_{ij}$. 因此, 内部求和 $\sum_{i=1}^n a_{ij} q_i$ 恰好是向量 $A^* Q$ 的第 $j$ 个分量. 

$$\langle L_2(P), Q\rangle = \sum_{j=1}^n \langle p_j, (A^* Q)_j\rangle = \langle P, A^* Q\rangle$$

### 求解 $L_1(P) = DP(x)Ax$ 的共轭
向量 $DP(x)Ax$ 的第 $i$ 个分量由下式给出:

$$(L_1(P))_i = \sum_{j=1}^n \frac{\partial p_i}{\partial x_j} (Ax)_j = \sum_{j=1}^n \sum_{k=1}^n a_{jk} x_k \frac{\partial p_i}{\partial x_j}$$

现在我们计算内积:

$$\langle L_1(P), Q\rangle = \sum_{i=1}^n \langle (L_1(P))_i, q_i\rangle = \sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n a_{jk} \langle x_k \frac{\partial p_i}{\partial x_j}, q_i\rangle$$

应用 **性质 A** 和 **性质 B**:
$\langle x_k \frac{\partial p_i}{\partial x_j}, q_i\rangle = \langle \frac{\partial p_i}{\partial x_j}, \frac{\partial q_i}{\partial x_k}\rangle = \langle p_i, x_j \frac{\partial q_i}{\partial x_k}\rangle$

将其代回求和式中:

$$\langle L_1(P), Q\rangle = \sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n a_{jk} \langle p_i, x_j \frac{\partial q_i}{\partial x_k}\rangle = \sum_{i=1}^n \langle p_i, \sum_{k=1}^n \frac{\partial q_i}{\partial x_k} (\sum_{j=1}^n a_{jk} x_j)\rangle$$

注意最内层求和 $\sum_{j=1}^n a_{jk} x_j = \sum_{j=1}^n (A^*)_{kj} x_j = (A^* x)_k$. 替换后得到:

$$\langle L_1(P), Q\rangle = \sum_{i=1}^n \langle p_i, \sum_{k=1}^n \frac{\partial q_i}{\partial x_k} (A^* x)_k\rangle$$

项 $\sum_{k=1}^n \frac{\partial q_i}{\partial x_k} (A^* x)_k$ 恰好是向量 $DQ(x)A^* x$ 的第 $i$ 个分量. 因此:

$$\langle L_1(P), Q\rangle = \sum_{i=1}^n \langle p_i, (DQ(x)A^* x)_i\rangle = \langle P, DQ(x)A^* x\rangle.$$

这完成了 $(ad_A)^* = ad_{A^*}$ 的证明.

## 映射: $(ad_A)^*=ad_{A^*}$ 的证明

只需证明: 算子 $L_A: P(y) \mapsto P(Ay)$ 的共轭为 $L_{A^*}: Q(y) \mapsto Q(A^*y)$

为了证明 $(L_A)^* = L_{A^*}$, 我们需要证明对于任意 $P, Q \in H_n^k$, 满足等式:
$$\langle P(Ay), Q(y)\rangle = \langle P(y), Q(A^*y)\rangle.$$

对于标量多项式 $f(y),g(y)$, 上面的等式退化为:
$$\langle f(Ay), g(y)\rangle = \langle f(y), g(A^*y)\rangle.$$

考虑等式右边. 引入中间变量 $u = A^*y$ (在实空间下 $A^* = A^T$). 故 $u_k = \sum_{j=1}^n a_{jk} y_j$.
根据多元复合函数的链式法则, 偏导数算子 $\nabla_y$ 作用于函数 $g(u)$ 时满足:
$$\frac{\partial}{\partial y_i} = \sum_{k=1}^n \frac{\partial u_k}{\partial y_i} \frac{\partial}{\partial u_k}=\sum_{k=1}^n a_{ik}\frac{\partial}{\partial u_k}$$
写成向量算子形式, 即:
$$\nabla_y = A \nabla_u$$

将上述算子关系代入右边内积定义:
$$\langle f(y), g(A^*y)\rangle = \left[ f(\nabla_y) g(A^*y) \right]_{y=0} = \left[ f(A \nabla_u) g(u) \right]_{u=0}.$$

以上证明了标量形式的结论. 自然地, 标量形式的结论的结论可推广到矢量形式.