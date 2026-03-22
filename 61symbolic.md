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