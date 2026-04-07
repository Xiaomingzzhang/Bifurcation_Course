# 数值方法

这一小节中, 我们介绍一些常见的数值方法. 尽管许多方法已经尽为人知了, 但对于初学者来说, 熟悉其原理与程序实现的细节也是十分关键的. 程序的实现仅取决于作者所熟悉的编程语言. 

## 低维系统的分岔曲线与分岔曲面
对于自治的光滑系统, 当系统的维数加上参数的个数小于等于 3, 我们可以在三维空间中将平衡点曲线或者平衡点曲面可视化. 绘制这些曲面或者曲线的函数是等值面或等值线函数, 如 Mathematica 中的 `ContourPlot` 和 `ContourPlot3D`. 另外, 由于平衡点的稳定性可以由偏导数的正负来判定, 我们可根据这个条件对曲线进行着色.
```mathematica
ClearAll[PlotBifurcation];
PlotBifurcation[f_, x_Symbol, mu_Symbol, muRange_List, xRange_List] :=
  Module[{dfdx, stablePlot, unstablePlot},(*1. 计算导数*)dfdx = D[f, x];
  (*2. 绘制稳定分支 (df/dx<0)*)
  stablePlot = 
   ContourPlot[
    f == 0, {mu, muRange[[1]], muRange[[2]]}, {x, xRange[[1]], 
     xRange[[2]]}, 
    ContourStyle -> {Blue, Thickness[0.005]},(*使用替换确保计算正确*)
    RegionFunction -> 
     Function[{mVal, xVal, z}, (dfdx /. {mu -> mVal, x -> xVal}) < 0],
     PlotPoints -> 50, MaxRecursion -> 3];
  (*3. 绘制不稳定分支 (df/dx>0)*)
  unstablePlot = 
   ContourPlot[
    f == 0, {mu, muRange[[1]], muRange[[2]]}, {x, xRange[[1]], 
     xRange[[2]]}, ContourStyle -> {Red, Dashed, Thickness[0.005]}, 
    RegionFunction -> 
     Function[{mVal, xVal, z}, (dfdx /. {mu -> mVal, x -> xVal}) > 0],
     PlotPoints -> 50, MaxRecursion -> 3];
  (*4. 合并图形并修正 PlotRange 优先级*)
  Show[{stablePlot, unstablePlot}, Frame -> True, 
   FrameLabel -> {ToString[mu], ToString[x]}, 
   PlotLabel -> "Bifurcation Diagram", GridLines -> Automatic, 
   GridLinesStyle -> 
    Directive[Gray, Dashed, Opacity[0.5]],
   PlotRange -> {muRange, xRange}, ImageSize -> Medium, 
   AspectRatio -> 0.8]]
Clear[x, mu];
GraphicsGrid[{{PlotBifurcation[\[Mu] + x^2, x, \[Mu], {-2, 2}, {-2, 2}],
    PlotBifurcation[\[Mu]*x - x^2, 
    x, \[Mu], {-2, 2}, {-2, 2}]}, {PlotBifurcation[\[Mu]*x - x^3, 
    x, \[Mu], {-2, 2}, {-2, 2}], 
   PlotBifurcation[\[Mu] + x - x^3, x, \[Mu], {-2, 2}, {-2, 2}]}}]
```
```{figure} ./asserts/figs/62_b_plot.png
:alt: 图片无法加载
:width: 100%
:align: center
```
上述代码的 matlab 实现 (由 AI 实现代码翻译):
```matlab
% 1. 环境准备
clear; clc;
syms mu x
muRange = [-2, 2];
xRange = [-2, 2];

% 2. 创建画布
figure('Name', '四种经典分岔图', 'Color', 'w', 'Position', [100, 100, 800, 600]);
tiledlayout(2, 2, 'TileSpacing', 'Compact');

% --- 案例 1: 鞍结分岔 (Saddle-node) ---
nexttile;
PlotBifurcation(mu + x^2, x, mu, muRange, xRange);
title('鞍结分岔: \mu + x^2 = 0');

% --- 案例 2: 跨临界分岔 (Transcritical) ---
nexttile;
PlotBifurcation(mu*x - x^2, x, mu, muRange, xRange);
title('跨临界分岔: \mu x - x^2 = 0');

% --- 案例 3: 音叉分岔 (Pitchfork) ---
nexttile;
PlotBifurcation(mu*x - x^3, x, mu, muRange, xRange);
title('音叉分岔: \mu x - x^3 = 0');

% --- 案例 4: 滞后分岔 (Imperfection) ---
nexttile;
PlotBifurcation(mu + x - x^3, x, mu, muRange, xRange);
title('滞后分岔: \mu + x - x^3 = 0');

function PlotBifurcation(f, x_sym, mu_sym, muRange, xRange)
    % 1. 计算偏导数 df/dx 用于判断稳定性
    dfdx_sym = diff(f, x_sym);

    % 将符号表达式转换为数值函数句柄，加快绘图速度
    % 'Vars', [mu_sym, x_sym] 确保输入顺序
    f_num = matlabFunction(f, 'Vars', [mu_sym, x_sym]);
    df_num = matlabFunction(dfdx_sym, 'Vars', [mu_sym, x_sym]);

    % 2. 定义包装函数：只有满足稳定性条件时才返回 f 的值，否则返回 NaN
    % 稳定分支: df/dx < 0
    f_stable = @(m, x) arrayfun(@(mv, xv) validate_branch(mv, xv, f_num, df_num, 'stable'), m, x);
    % 不稳定分支: df/dx > 0
    f_unstable = @(m, x) arrayfun(@(mv, xv) validate_branch(mv, xv, f_num, df_num, 'unstable'), m, x);

    % 3. 绘制图形
    hold on;
    
    % 绘制稳定分支 (蓝色实线)
    h1 = fimplicit(f_stable, [muRange xRange], 'Color', 'b', 'LineWidth', 2);
    
    % 绘制不稳定分支 (红色虚线)
    h2 = fimplicit(f_unstable, [muRange xRange], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

    % 4. 美化图形
    grid on;
    ax = gca;
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.5;
    xlabel(char(mu_sym));
    ylabel(char(x_sym));
    
    % 绘制坐标轴辅助线
    line(muRange, [0 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    line([0 0], xRange, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    
    hold off;
end

% 辅助函数：根据导数正负过滤分支
function val = validate_branch(mu_val, x_val, f_func, df_func, type)
    try
        d = df_func(mu_val, x_val);
        if strcmp(type, 'stable') && d < 0
            val = f_func(mu_val, x_val);
        elseif strcmp(type, 'unstable') && d > 0
            val = f_func(mu_val, x_val);
        else
            val = NaN;
        end
    catch
        % 处理导数为常数或计算异常的情况
        val = f_func(mu_val, x_val); 
    end
end
```

尖点分岔的分岔曲面:
```mathematica
ContourPlot3D[
 c + d x - x^3 == 0, {c, -20, 20}, {d, -20, 20}, {x, -10, 10}, 
 MeshStyle -> None, PlotPoints -> 200, AxesLabel -> {"c", "d", "x"}, 
 ColorFunction -> Function[{c, d, x, f}, If[d - 3 x^2 >= 0, Red, Blue]],
  ColorFunctionScaling -> False, 
 ContourStyle -> Directive[Opacity[0.8], Specularity[White, 30]]]
```
```{figure} ./asserts/figs/62_b_plot_3d.png
:alt: 图片无法加载
:width: 80%
:align: center
```
上述代码的 matlab 实现 (由 AI 实现代码翻译):
```matlab
% 定义范围
cRange = [-20, 20];
dRange = [-20, 20];
xRange = [-10, 10];

figure('Color', 'w');
hold on;

% 1. 绘制稳定分支 (d - 3*x^2 < 0)
% 我们通过逻辑判断将不符合条件的区域设为 NaN
f_stable = @(c, d, x) deal_with_condition(c, d, x, @(c,d,x) d - 3*x.^2 < 0);
h_stable = fimplicit3(f_stable, [cRange dRange xRange], ...
    'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.8,'MeshDensity', 100);

% 2. 绘制不稳定分支 (d - 3*x^2 >= 0)
f_unstable = @(c, d, x) deal_with_condition(c, d, x, @(c,d,x) d - 3*x.^2 >= 0);
h_unstable = fimplicit3(f_unstable, [cRange dRange xRange], ...
    'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.8,'MeshDensity', 100);

% 3. 图形美化
xlabel('c'); ylabel('d'); zlabel('x');
title('尖点灾变曲面 (Cusp Catastrophe)');
view(35, 30);
grid on;

% 设置光照效果，模拟 Mathematica 的 Specularity
camlight headlight;
lighting gouraud;
material shiny;

hold off;

% 辅助函数：根据稳定性条件过滤曲面
function val = deal_with_condition(c, d, x, condition_func)
    % 隐函数方程：c + d*x - x^3 = 0
    f_val = c + d.*x - x.^3;
    
    % 检查条件
    mask = condition_func(c, d, x);
    
    % 不满足条件的地方返回 NaN，fimplicit3 将不会绘制这些点
    f_val(~mask) = NaN;
    val = f_val;
end
```
## 单自由度保守系统的相图


## 使用自动微分数值计算多重线性型

## 变分方程与牛顿迭代法计算周期轨道
