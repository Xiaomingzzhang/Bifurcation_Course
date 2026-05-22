# 符号计算

我们将介绍使用 Mathematica 软件, 计算分岔中涉及到的繁杂的计算, 如规范型的系数, 中心流形的系数, 投影法的应用等.
## Hopf 分岔 $c_{1}(\mu)$ 的计算
根据 Hopf 分岔的范式推导, 其分岔主要由系数 $c_{1}(\mu)$ 来确定, 我们介绍如何计算.
### 方法一: 对逆变换求导
第一种方法是对逆变换
$$
\label{61inverse_transform}
w=z-h_{2}(z,\bar{z},\mu)+O(|z|^3)
$$
求导. 然后代入
$$
\dot{z}=\lambda z+g_{2}(z,\bar{z},\mu)+g_{3}(z,\bar{z},\mu)+O(|z|^4),
$$
再将式子中的 $z$ 用
$$
z=w+h_{2}(w,\bar{w},\mu)
$$
替换. 然而, 我们发现, 最终得到的三次项系数里面包含逆变换 [](#61inverse_transform) 的三次项的系数. 如果运用这种策略, 我们首先应该得到逆变换 [](#61inverse_transform) 的三次项.

首先我们定义几个函数:
```mathematica
ClearAll["Global`*"];
GeneratePoly[w_, wb_, g_, l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,
   0}}*)indices = FrobeniusSolve[{1, 1}, l];
   (*2. 构造每一项：(1/(j! k!))*h[j,k]*w^j*wb^k*)
   terms = Table[
     With[{j = idx[[1]], k = idx[[2]]}, (1/(j!  k!))*
       Subscript[g, j, k]*w^j*wb^k], {idx, indices}];
   (*3. 求和并返回*)Total[terms]];
GenerateRules[l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,
   0}}*)indices = FrobeniusSolve[{1, 1}, l];
   (*2. 替换规则*)
   terms = Table[
     With[{j = idx[[1]], 
       k = idx[[2]]}, {Subscript[h, j, k] -> Subscript[g, j, k]/(
        j \[Lambda] + k \[Lambda]b - \[Lambda]), 
       Subscript[hb, j, k] -> Subscript[gb, j, k]/(
        j \[Lambda]b + k \[Lambda] - \[Lambda]b)}], {idx, indices}];
   (*3. 返回替换规则*)terms // Flatten];
SequentialReplace[expr_, rules_List] := Fold[ReplaceAll, expr, rules];
```
由于符号计算软件对于符号的共轭运算很不智能, 我们只能手动处理共轭情况. 代码中**若 $a$ 是变量, 那么变量 $ab$ 就表示其共轭.**

上面的函数 `GeneratePoly` 是为了生成 $h_{2}(z,\bar{z},\mu)$ 这样的多项式. 例如, 输入 `GeneratePoly[z,zb,h,2]` 我们得到
$$
\frac{1}{2} z^2 h_{2,0}+z \text{zb} h_{1,1}+\frac{1}{2} \text{zb}^2 h_{0,2},
$$
输入 `GeneratePoly[z,zb,g,4]` 我们得到
$$
\frac{1}{24} z^4 g_{4,0}+\frac{1}{6} z^3 \text{zb} g_{3,1}+\frac{1}{4} z^2 \text{zb}^2 g_{2,2}+\frac{1}{6} z \text{zb}^3 g_{1,3}+\frac{1}{24} \text{zb}^4 g_{0,4}.
$$
`GenerateRules` 是按照推导得到的替换规则, `SequentialReplace` 可以对表达式施加连续的替换.

首先我们利用展开的方式求解逆变换的三次项系数, 三次项的系数为 $p$:
```mathematica
(*计算逆变换w=z-h2(z,zb)+p3(z,zb)的三次项的系数*)
(*三次项系数*)
indices = FrobeniusSolve[{1, 1}, 3];
(*z=w+h2(w,wb)中代入w=z-h2(z,zb)+p3(z,zb)*)
expression = 
  z - w - GeneratePoly[w, wb, h, 2] /. {w -> 
      z - GeneratePoly[z, zb, h, 2] + GeneratePoly[z, zb, p, 3], 
     wb -> zb - GeneratePoly[zb, z, hb, 2] + 
       GeneratePoly[zb, z, pb, 3]} // Expand;
(*求解满足i+j=3的p[i,j]的系数*)
solvelist = 
  Table[With[{i = idx[[1]], j = idx[[2]]}, Subscript[p, i, j]], {idx, indices}];
coep = Table[
   With[{i = idx[[1]], j = idx[[2]]}, 
    Coefficient[expression, z^i zb^j] /. {z -> 0, zb -> 0}], {idx, 
    indices}];
sol = Solve[# == 0 & /@ coep, solvelist][[1]].
```
我们得到
$$
\left\{p_{0,3}\to 3 \left(h_{0,2} \text{hb}_{2,0}+h_{0,2} h_{1,1}\right),p_{1,2}\to h_{1,1} \text{hb}_{2,0}+2 h_{0,2} \text{hb}_{1,1}+2 h_{1,1}^2+h_{0,2} h_{2,0},p_{2,1}\to h_{0,2} \text{hb}_{0,2}+2 h_{1,1} \text{hb}_{1,1}+3 h_{1,1} h_{2,0},p_{3,0}\to 3 \left(h_{1,1} \text{hb}_{0,2}+h_{2,0}^2\right)\right\}
$$
现在我们可以采用刚刚讲述过的策略了, 对逆变换求导, 而后再逐步替换:
```mathematica
(*z,zb的替换规则*)
rule1 = z -> (w + GeneratePoly[w, wb, h, 2]);
rule2 = zb -> (wb + GeneratePoly[wb, w, hb, 2]);
(*z的导数,zb的导数的替换规则*)
rule3 = dz -> (\[Lambda]  z + GeneratePoly[z, zb, g, 2] + 
     GeneratePoly[z, zb, g, 3]);
rule4 = cojdz -> (\[Lambda]b  zb + GeneratePoly[zb, z, gb, 2] + 
     GeneratePoly[zb, z, gb, 3]);
preresult = (dz - D[GeneratePoly[z, zb, h, 2], z]*dz - 
    D[GeneratePoly[z, zb, h, 2], zb]*cojdz + 
    D[GeneratePoly[z, zb, p, 3], z]*dz + 
    D[GeneratePoly[z, zb, p, 3], zb]*cojdz);

(*对w=z-h2(z,zb)+p3(z,zb)求导后逐步替换*)
result = 
  SequentialReplace[
    preresult, {sol, rule4, rule3, rule2, rule1, GenerateRules[2]}] //
    Expand;
(*验证二次项系数,应该全部为0,提取三次项系数*)
indices2 = FrobeniusSolve[{1, 1}, 2];
Table[With[{i = idx[[1]], j = idx[[2]]}, 
   Coefficient[result, w^i wb^j] /. {w -> 0, wb -> 0}], {idx, 
   indices2}] // Simplify
Coefficient[result, w^2 wb]
```
运算结果应该得到二次项系数的三个 0, 以及所需要的三次项 `w^2 wb` 的系数:
$$
\frac{g_{0,2} \text{gb}_{0,2}}{2 (2 \lambda -\text{$\lambda $b})}+\frac{\lambda  g_{0,2} \text{gb}_{0,2}}{(2 \lambda -\text{$\lambda $b}) (2 \text{$\lambda $b}-\lambda )}-\frac{g_{0,2} \text{gb}_{0,2}}{2 (2 \text{$\lambda $b}-\lambda )}-\frac{\text{$\lambda $b} g_{0,2} \text{gb}_{0,2}}{2 (2 \lambda -\text{$\lambda $b}) (2 \text{$\lambda $b}-\lambda )}+\frac{g_{1,1} \text{gb}_{1,1}}{\lambda }+\frac{g_{1,1} g_{2,0}}{2 \lambda }+\frac{g_{1,1} g_{2,0}}{\text{$\lambda $b}}+\frac{g_{2,1}}{2}.
$$
中间的几项:
$$
\frac{\lambda  g_{0,2} \text{gb}_{0,2}}{(2 \lambda -\text{$\lambda $b}) (2 \text{$\lambda $b}-\lambda )}-\frac{g_{0,2} \text{gb}_{0,2}}{2 (2 \text{$\lambda $b}-\lambda )}-\frac{\text{$\lambda $b} g_{0,2} \text{gb}_{0,2}}{2 (2 \lambda -\text{$\lambda $b}) (2 \text{$\lambda $b}-\lambda )}
$$
实际为 $0$. 那么最终结果即:
$$
\frac{g_{0,2} \text{gb}_{0,2}}{2 (2 \lambda -\text{$\lambda $b})}+\frac{g_{1,1} \text{gb}_{1,1}}{\lambda }+\frac{g_{1,1} g_{2,0}}{2 \lambda }+\frac{g_{1,1} g_{2,0}}{\text{$\lambda $b}}+\frac{g_{2,1}}{2}.
$$
### 方法二: 对原变换求导

上面的方法实际上需要求逆变换的三次项系数, 较为麻烦. 实际上, 我们可以对原变换
$$
z=w+h_{2}(w,\bar{w},\mu)
$$
求导. 求完导以后, 方程的左边用
$$
\dot{z}=\lambda z+g_{2}(z,\bar{z},\mu)+g_{3}(z,\bar{z},\mu)+O(|z|^4),
$$
替换, 方程的右边中关于 $w$ 的导数用
$$
\dot{w}=\lambda w+c_{1}w^2\bar{w}+(|z|^4)
$$
替换. 自然地, 上面的关于 $w$ 的微分方程有如此简单的形式是因为我们已经取定了 $h_{2}$ 的系数使得全部的二阶项可以消去. 最后再将方程的左边关于 $z$ 的项用原变换替换即可. 这个方法比上面的方法要简单许多, 完整的代码如下:
```mathematica
ClearAll["Global`*"];
GeneratePoly[w_, wb_, g_, l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,
   0}}*)indices = FrobeniusSolve[{1, 1}, l];
   (*2. 构造每一项：(1/(j! k!))*h[j,k]*w^j*wb^k*)
   terms = Table[
     With[{j = idx[[1]], k = idx[[2]]}, (1/(j!  k!))*
       Subscript[g, j, k]*w^j*wb^k], {idx, indices}];
   (*3. 求和并返回*)Total[terms]];
GenerateRules[l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,
   0}}*)indices = FrobeniusSolve[{1, 1}, l];
   (*2. 替换规则*)
   terms = Table[
     With[{j = idx[[1]], 
       k = idx[[2]]}, {Subscript[h, j, k] -> Subscript[g, j, k]/(
        j \[Lambda] + k \[Lambda]b - \[Lambda]), 
       Subscript[hb, j, k] -> Subscript[gb, j, k]/(
        j \[Lambda]b + k \[Lambda] - \[Lambda]b)}], {idx, indices}];
   (*3. 返回替换规则*)terms // Flatten];
SequentialReplace[expr_, rules_List] := 
  Fold[ReplaceAll, expr, rules];
(*z,zb,z的导数的替换规则*)
rule1 = z -> (w + GeneratePoly[w, wb, h, 2]);
rule2 = zb -> (wb + GeneratePoly[wb, w, hb, 2]);
rule3 = dz -> (\[Lambda]  z + GeneratePoly[z, zb, g, 2] + 
     GeneratePoly[z, zb, g, 3]);
rule4 = dw -> \[Lambda]  w + c  w^2  wb;
rule5 = dwb -> \[Lambda]b  wb + cb  w  wb^2;
preresult = (dz - dw - D[GeneratePoly[w, wb, h, 2], w]*dw - 
    D[GeneratePoly[w, wb, h, 2], wb]*dwb);
(*对z=w+h2(w,wb)求导后逐步替换*)
result = 
  SequentialReplace[
    preresult, {rule4, rule5, rule3, rule2, rule1, 
     GenerateRules[2]}] // Expand;
(*验证二次项系数,应该全部为0,提取三次项系数*)
indices2 = FrobeniusSolve[{1, 1}, 2];
Table[With[{i = idx[[1]], j = idx[[2]]}, 
   Coefficient[result, w^i wb^j] /. {w -> 0, wb -> 0}], {idx, 
   indices2}] // Simplify
Coefficient[result, w^2 wb] /. {w -> 0, wb -> 0}
```
运算得到二次项系数全为 $0$, 并得到 `w^2 wb` 的系数为:
$$
\frac{g_{0,2} \text{gb}_{0,2}}{2 (2 \lambda -\text{$\lambda $b})}+\frac{g_{1,1} \text{gb}_{1,1}}{\lambda }+\frac{g_{1,1} g_{2,0}}{2 \lambda }+\frac{g_{1,1} g_{2,0}}{\text{$\lambda $b}}+\frac{g_{2,1}}{2}-c.
$$

以上代码的计算思路可以进行推广, 应用到其他分岔的范式系数的计算.

## NS 分岔 $c_{1}(\mu)$ 的计算

与之前的思路类似, 这里我们给出代码:
```mathematica
ClearAll["Global`*"];
GeneratePoly[w_, wb_, g_, l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,0}}*)
   indices = FrobeniusSolve[{1, 1}, l];
   (*2. 构造每一项：(1/(j! k!))*h[j,k]*w^j*wb^k*)
   terms = Table[
     With[{j = idx[[1]], k = idx[[2]]}, (1/(j! k!))*Subscript[g, j, k]*
       w^j*wb^k], {idx, indices}];
   (*3. 求和并返回*)Total[terms]];
GenerateRules[l_Integer] := 
  Module[{indices, terms},(*1. 找到所有满足 j+k=
   l 的非负整数解*)(*FrobeniusSolve[{1,1},l] 会直接返回 {{0,l},{1,l-1},...,{l,0}}*)
   indices = FrobeniusSolve[{1, 1}, l];
   (*2. 替换规则*)
   terms = Table[
     With[{j = idx[[1]], 
       k = idx[[2]]}, {Subscript[h, j, k] -> 
        Subscript[g, j, k]/(\[Lambda] - \[Lambda]b^k*\[Lambda]^j), 
       Subscript[hb, j, k] -> 
        Subscript[gb, j, 
         k]/(\[Lambda]b - \[Lambda]b^j*\[Lambda]^k)}], {idx, indices}];
   (*3. 返回替换规则*)terms // Flatten];
SequentialReplace[expr_, rules_List] := Fold[ReplaceAll, expr, rules];
(*变换为:w=z+h2(z,zb), 计算逆变换z=w-h2(w,wb)+p3(w,wb)的三次项的系数*)
(*三次项系数*)
indices = FrobeniusSolve[{1, 1}, 3];
(*w=z+h2(z,zb)中代入z=w-h2(w,wb)+p3(w,wb)*)
expression = 
  w - z - GeneratePoly[z, zb, h, 2] /. {z -> 
      w - GeneratePoly[w, wb, h, 2] + GeneratePoly[w, wb, p, 3], 
     zb -> wb - GeneratePoly[wb, w, hb, 2] + 
       GeneratePoly[wb, w, pb, 3]} // Expand;
(*求解满足i+j=3的p[i,j]的系数*)
solvelist = 
  Table[With[{i = idx[[1]], j = idx[[2]]}, Subscript[p, i, j]], {idx, 
    indices}];
coep = Table[
   With[{i = idx[[1]], j = idx[[2]]}, 
    Coefficient[expression, w^i wb^j] /. {w -> 0, wb -> 0}], {idx, 
    indices}];
sol = Solve[# == 0 & /@ coep, solvelist][[1]];
solb = sol /. {p -> pb, h -> hb, hb -> h};
rule1 = z -> (w - GeneratePoly[w, wb, h, 2] + 
     GeneratePoly[w, wb, p, 3]);
rule2 = zb -> (wb - GeneratePoly[wb, w, hb, 2] + 
     GeneratePoly[wb, w, pb, 3]);
rule3 = zprime -> (\[Lambda] z + GeneratePoly[z, zb, g, 2] + 
     GeneratePoly[z, zb, g, 3]);
rule4 = zbprime -> (\[Lambda]b zb + GeneratePoly[zb, z, gb, 2] + 
     GeneratePoly[zb, z, gb, 3]);
preresult = zprime + GeneratePoly[zprime, zbprime, h, 2];
result = 
  SequentialReplace[
   preresult, {rule3, rule4, rule2, rule1, sol, solb, GenerateRules[2]}];
(*验证二次项系数,应该全部为0,提取三次项系数*)
indices2 = FrobeniusSolve[{1, 1}, 2];
Table[With[{i = idx[[1]], j = idx[[2]]}, 
   Coefficient[result, w^i wb^j] /. {w -> 0, wb -> 0}], {idx, 
   indices2}] // Simplify
Coefficient[result, w^2 wb] /. {w -> 0, wb -> 0} // Simplify
```
得到二次项系数全部为 0, 三次项 $w^2\bar{w}$ 的系数为:
$$
\frac{1}{2} \left(\frac{\text{$\lambda $b} g_{0,2} \text{gb}_{0,2}}{\lambda -\text{$\lambda $b}^2}-\frac{g_{1,1} \left(g_{2,0} (3 \lambda  \text{$\lambda $b}-2 \lambda -\text{$\lambda $b})+2 (\lambda -1) \lambda  \text{gb}_{1,1}\right)}{(\lambda -1) \lambda  (\text{$\lambda $b}-1)}+g_{2,1}\right)
$$

```{prf:remark}
有些教材上 $c_1(\mu)$ 的表达式与上式不同, 实际上, 由于 $c_1(\mu)$ 仅在分岔点处的计算有意义, 利用 $|\lambda|=1$ 即可得到等价的式子.
```

## 多重线性型的计算
在投影法中, 经常要计算矢量场的多重线性型, 特别是当维数较高且矢量场表达式复杂时, 手工计算非常繁琐且很容易出错.

回顾由多重线性型所定义的多维多元函数的泰勒展开. 如果 $\mathbf{f}:\mathbb{R}^{n}\rightarrow\mathbb{R}^m$, $\mathbf{f}$ 的形式如下:
$$
\mathbf{f}(\mathbf{x})=\Big(f_1(\mathbf{x}),f_2(\mathbf{x}),\cdots,f_m(\mathbf{x}) \Big)^{T}.
$$
 $\mathbf{f}:\mathbb{R}^{n}\rightarrow\mathbb{R}^m$ 是充分光滑的, 那么在 $\mathbf{a}\in \mathbb{R}^{n}$ 处,
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

我们要设计一个程序函数, 这个函数的调用变量为:
- 一个 $m$ 维矢量场;
- 矢量场的状态变量;
- 泰勒展开的出发点;
- 所求的多重线性型的阶次 $l$.

这个函数的输出为一个阶次为 $l$ 的多重线性型, 其也为一个函数, 接收 $l$ 个 $n$ 维矢量, 返回一个 $m$ 维矢量.

使用 Mathematica, 我们仅需要几行代码便可以实现上面的复杂逻辑.
```mathematica
ClearAll["Global`*"];
(*f:向量场或映射 (List) x:自变量列表 (List) a:展开点 (List) k:多重线性型的阶数 (Integer)*)
MultilinearForm[f_, x_, a_, k_] := 
  Module[{tensorD}, tensorD = D[f, {x, k}] /. Thread[x -> a];
   Function[Null, 
    Block[{vecs = {##}}, Fold[Dot, tensorD, Reverse[vecs]]]]];
```
上面函数的代码逻辑如下. 首先, 对于一个以 $n$ 个自变量的 $m$ 维的向量函数
$$
\mathbf{f}(\mathbf{x})=\Big(f_1(\mathbf{x}),f_2(\mathbf{x}),\cdots,f_m(\mathbf{x}) \Big)^{T},
$$
`D[f[x],{x,l}]` 生成一个
$$
\{m,\underbrace{n,n,\cdots,n}_{l\text{个}} \}
$$ 维的张量. 具体地, 如果对这个张量进行指标索引 $k,i_1,\cdots,i_{l}$, 得到的是数:
$$
\frac{\partial ^lf_k}{\partial x_{i_1}\partial x_{i_2}\cdots\partial x_{i_l}}(\mathbf{a}).
$$
这正是我们想要的 $l$ 阶多重线性型的系数. 两个张量 $T,Q$ 进行缩并 `Dot` 运算也就是 $T$ 的最后一个指标与 $Q$ 的第一个指标联合起来求和. 因此如果张量 `D[f[x],{x,l}]` 与矢量 $\xi$ 进行 `Dot` 运算, 那么就得到一个阶数为
$$
\{m,\underbrace{n,n,\cdots,n}_{l-1\text{个}} \}
$$
的张量, 它的第 $k,i_1,\cdots,i_{l-1}$ 个变量为:
$$
\sum_{i_{l}=1}^{n}\frac{\partial ^lf_k}{\partial x_{i_1}\partial x_{i_2}\cdots\partial x_{i_l}}(\mathbf{a})\xi_{i_l}.
$$

所以, 对于 $l$ 个矢量 `v_{1},v_{2},...,v_{l}`, 如果我们将 `v_{l}` 先与 `D[f[x],{x,l}]` 进行缩并, 得到的结果再与 `v_{l-1}` 缩并, 一直重复下去直到与 `v_{1}` 缩并得到一个 $m$ 维矢量, 即得到这个 $l$ 阶的多重线性型对于 `v_{1},v_{2},...,v_{l}` 的作用结果.

对于课本上的干摩擦的例子 (p97):
$$
\dot{x}&=y,\\
\dot{y}&=-\frac{k}{m}y+F(y-v)/m.
$$
运行代码
```mathematica
vars = {x, y};
f = {y, -\[Omega]^2 x + F[y - v]/m};
a = {F[-v]/m, 0}; 
A = MultilinearForm[f, vars, a, 1];
B = MultilinearForm[f, vars, a, 2];
G = MultilinearForm[f, vars, a, 3];
var\[Xi] = {Subscript[\[Xi], 1], Subscript[\[Xi], 2]};
var\[Eta] = {Subscript[\[Eta], 1], Subscript[\[Eta], 2]};
var\[Zeta] = {Subscript[\[Zeta], 1], Subscript[\[Zeta], 2]};
A[var\[Xi]]
B[var\[Xi], var\[Eta]]
G[var\[Xi], var\[Eta], var\[Zeta]]
```
得到
$$
\left\{\xi _2,\frac{\xi _2 F'(-v)}{m}-\xi _1 \omega ^2\right\},\\
\left\{0,\frac{\eta _2 \xi _2 F''(-v)}{m}\right\},\\
\left\{0,\frac{\zeta _2 \eta _2 \xi _2 F^{(3)}(-v)}{m}\right\},
$$
与书本中公式相同.

## 中心流形的计算: 坐标变换法

下面我们给出一个例子, 来看如何使用程序计算中心流形.

考虑 Lorenz 系统:
$$
\begin{aligned}
\dot{x}& = \sigma (y - x), \\
\dot{y}& = x(\rho - z) - y, \\
\dot{z}& = xy - \beta z.
\end{aligned}
$$
这里将取 $\rho=1$. 首先我们定义矢量场, 并得到系统在原点的 Jacobi 矩阵的特征值与特征矢量:
```mathematica
f[x_, y_, z_] := {\[Sigma] (y - x), x - y - x z, -\[Beta] z + x y};
A = Grad[f[x, y, z], {x, y, z}] /. {x -> 0, y -> 0, z -> 0};
Eigensystem[A]
```
我们得到:
$$
\left(
\begin{array}{ccc}
 0 & -\beta  & -\sigma -1 \\
 \{1,1,0\} & \{0,0,1\} & \{-\sigma ,1,0\} \\
\end{array}
\right)
$$
这样特征值相应的特征矢量排成一列就是我们要的线性变换矩阵, 这样定义新系统为:
```mathematica
P = Transpose[Eigensystem[A][[2]]];
newf[u_, v_, w_] = (Inverse[P] . f[x, y, z]) /. 
   Thread[{x, y, z} -> P . {u, v, w}];
```

下面我们假设中心流形的形式, 并代入到不变方程中. 要计算到直到 $k$ 阶的中心流形, 其一般形式是一个 $n$ 维的 $m$ 个变量的直到 $k$ 阶的 $n$ 维多项式. 对于我们这个例子, $n=2,m=1$, 并取 $k=5$. 为此我们定义一个一般的函数, 能生成任意的这种多项式, 这样程序稍作修改即可计算任意形式的中心流形:
```mathematica
(*定义一个函数,它生成n维m个变量的,从 k1 阶到 k2 阶的多项式向量*)
h[n_, m_, k1_, k2_, var_] := 
  Module[{indices, terms, polys, coeffs = {}},(*1. 找到所有满足 k1<=d_1+...+
   d_m<=k2 的非负整数解*)(*遍历所有的目标阶数 s (从 k1 到 k2),对每个 s 找到精确解并展平合并*)
   indices = 
    Flatten[Table[FrobeniusSolve[ConstantArray[1, m], s], {s, k1, k2}],
      1];
   polys = 
    Table[(*2. Construct each term for the i-th dimension*)
     terms = Table[
       With[{idx = currIdx},(*将当前生成的系数收集到 coeffs 列表中*)
        AppendTo[coeffs, Subscript[c, i, Sequence @@ idx]];
        (*Subscript[c,i,idx1,idx2...]*var[1]^idx1*var[2]^idx2...*)
        Subscript[c, i, Sequence @@ idx]*
         Product[var[[p]]^idx[[p]], {p, 1, m}]], {currIdx, indices}];
     (*3. Sum and return*)Total[terms], {i, 1, n}];
   {polys, coeffs}];
(*使得输入n,m,k,var时中心流形默认从2阶开始到k阶*)
h[n_, m_, k_, var_] := h[n, m, 2, k, var];
```
输入
```mathematica
h[2, 2, 3, {x, y}][[1]]
```
我们得到:
$$
\left\{x^3 c_{1,3,0}+x^2 y c_{1,2,1}+x^2 c_{1,2,0}+x y^2 c_{1,1,2}+x y c_{1,1,1}+y^3 c_{1,0,3}+y^2 c_{1,0,2},x^3 c_{2,3,0}+x^2 y c_{2,2,1}+x^2 c_{2,2,0}+x y^2 c_{2,1,2}+x y c_{2,1,1}+y^3 c_{2,0,3}+y^2 c_{2,0,2}\right\}
$$
这是一个 $2$ 维的以 $x,y$ 为变量的二阶到三阶的多维多项式, 且有待定系数, 这正是我们想要的函数.

再定义一个用于提取一个多项式系数, 并且使得多项式系数全为 $0$ 来求解待定系数的函数:
```mathematica
(*通用的系数提取与求解函数*)
SolveManifoldCoefficients[eq_, vars_, coeffs_, k_] := 
  Module[{rules, eqList, 
    sol},(*利用 CoefficientRules 提取每一维多项式中所有项的系数规则*)(*Expand 用于确保多项式被完全展\
开,避免提取遗漏*)rules = Flatten[CoefficientRules[Expand[#], vars] & /@ eq];
   (*过滤出所有总阶数小于等于 k 的项,并提取出其对应的系数表达式*)
   eqList = Cases[rules, ({n_} -> coef_) /; n <= k :> coef];
   (*令所有提取出的系数表达式等于 0,并对未知系数 coeffs 进行求解*)
   sol = Solve[Thread[eqList == 0], coeffs];
   sol[[1]]];
```

利用这两个函数我们就可以进行最终的中心流形的计算了:
```mathematica
order = 5;
centralmanifold = h[2, 1, order, {u}][[1]];
coeffs = h[2, 1, order, {u}][[2]];
eq = D[centralmanifold, 
       u]*(newf[u, v, w][[1]] /. Thread[{v, w} -> centralmanifold]) - 
     newf[u, v, w][[2 ;; 3]] /. Thread[{v, w} -> centralmanifold] // 
   Simplify;
sol = SolveManifoldCoefficients[eq, {u}, coeffs, order]
finalcentralmanifold = centralmanifold /. sol
reducedsystem = 
 Series[newf[u, v, w][[1]] /. 
   Thread[{v, w} -> finalcentralmanifold], {u, 0, order}]
```
得到一系列结果:
$$
\left\{c_{1,2}\to \frac{1}{\beta },c_{1,3}\to 0,c_{1,4}\to \frac{\beta  \sigma -\beta +2 \sigma ^2+2 \sigma }{\beta ^3 (\sigma +1)^2},c_{1,5}\to 0,c_{2,2}\to 0,c_{2,3}\to -\frac{1}{\beta  (\sigma +1)^2},c_{2,4}\to 0,c_{2,5}\to \frac{-5 \beta  \sigma +\beta -2 \sigma ^2-2 \sigma }{\beta ^3 (\sigma +1)^4}\right\}\\
\left\{\frac{u^4 \left(\beta  \sigma -\beta +2 \sigma ^2+2 \sigma \right)}{\beta ^3 (\sigma +1)^2}+\frac{u^2}{\beta },\frac{u^5 \left(-5 \beta  \sigma +\beta -2 \sigma ^2-2 \sigma \right)}{\beta ^3 (\sigma +1)^4}-\frac{u^3}{\beta  (\sigma +1)^2}\right\}\\
-\frac{\sigma  u^3}{\beta  (\sigma +1)}-\frac{\sigma  u^5 \left(2 \beta  \sigma -\beta +2 \sigma ^2+2 \sigma \right)}{\beta ^3 (\sigma +1)^3}+O\left(u^6\right)
$$
分别为中心流形的系数列表, 中心流形的表达式以及约化系统的矢量场.

由于函数 `h` 和 `SolveManifoldCoefficients` 的一般性, 上面的程序稍作修改可以求解任意有限维系统与任意阶次的中心流形, 前提是线性坐标变换可以得到 (有时求解特征值要求解高阶代数方程, 往往得不到线性坐标变换的解析表达式).

