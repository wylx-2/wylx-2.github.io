---
title: 《Riemann Solvers and Numerical Methods for Fluid Dynamics》读书笔记
author: wylx2
date: 2026-06-30 12:00:00 +0800
categories: [读书笔记, 流体力学]
tags: [读书笔记, CFD, Riemann, 流体力学, 双曲系统, Euler方程]
description: 关于 Riemann Solvers and Numerical Methods for Fluid Dynamics 一书的读书笔记，涵盖双曲型偏微分方程基本性质、Euler 方程特性及 Riemann 问题精确解法。
math: true
media_subpath: /assets/img/riemann-solvers/
---

本文是 E.F. Toro 所著《Riemann Solvers and Numerical Methods for Fluid Dynamics》的读书笔记，涵盖第 2–4 章的核心内容：双曲型 PDE 的基本理论、Euler 方程的性质、以及 Riemann 问题的精确解法。

## 1. 双曲型偏微分方程基本概念

### 1.1 拟线性系统

一阶偏微分系统可写作：

$$
\begin{equation}
\mathbf U_t +\mathbf {AU}_x +\mathbf B = \mathbf 0
\end{equation}
$$

- 若 $$\mathbf A$$ 和 $$\mathbf B$$ 为常数，则为**常系数线性系统**
- 若 $$\mathbf A$$ 和 $$\mathbf B$$ 仅依赖于 $$(x,t)$$，则为**变系数线性系统**
- 若 $$\mathbf B$$ 线性依赖于 $$\mathbf U$$，仍是线性系统；$$\mathbf B=\mathbf 0$$ 则称系统具有**齐次性**
- 若 $$\mathbf A = \mathbf A(\mathbf U)$$，则为**拟线性系统**（非线性系统）

**守恒律**形式：

$$
\mathbf U_t+\mathbf F(\mathbf U)_x=0
$$

Jacobi 矩阵：

$$
\mathbf{A(U)}=\partial \mathbf F/\partial \mathbf U=
\begin{bmatrix}
\partial f_1/\partial u_1 & \dots & \partial f_1/\partial u_m \\
\vdots & \vdots & \vdots \\
\partial f_m/\partial u_1 & \dots & \partial f_m/\partial u_m
\end{bmatrix}
$$

**特征值与特征向量**：

$$
\mathbf {AK}^{(i)}=\lambda_i\mathbf K^{(i)},\qquad
\mathbf L^{(i)}\mathbf A=\lambda_i\mathbf L^{(i)}
$$

**双曲性定义**：系统在 $$(x,t)$$ 处具有双曲性，当 $$\mathbf A$$ 具有 $$m$$ 个实数特征值和 $$m$$ 个线性无关的特征向量。特征值互异时称具有**严格双曲性**。（特征值均不为实数则为椭圆系统）

### 1.2 线性对流方程

考虑线性对流方程：

$$
u_t+au_x=0,\quad u(x,0)=u_0(x),\quad -\infty <x<\infty,\;t>0
$$

沿特征线 $$\text dx/\text dt=a$$，有 $$\text du/\text dt=0$$，即 $$u$$ 沿特征线不变。通解为：

$$
u(x,t)=u_0(x-at)
$$

波可视为以有限速度传播的扰动的可识别特征。

**Riemann 问题**：

$$
u(x,0)=u_0(x)=
\begin{cases}
u_L  & \text{if } x<0 \\
u_R  & \text{if } x>0
\end{cases}
$$

间断跟随特征线以波速 $$a$$ 移动。

### 1.3 线性双曲系统

考察：

$$
\mathbf U_t+\mathbf {AU}_x = \mathbf 0
$$

对角化 $$\mathbf A = \mathbf {K\Lambda K}^{-1}$$，引入特征变量 $$\mathbf {W} = \mathbf K^{-1}\mathbf U$$，解耦得到：

$$
\frac{\partial w_i}{\partial t}+\lambda_i\frac{\partial w_i}{\partial x}=0
$$

系统解为 $$m$$ 个波的叠加：

$$
\mathbf U(x,t)=\sum_{i=1}^{m}{w_i^{(0)}(x-\lambda_it)\mathbf K^{(i)}}
$$

**Riemann 问题的解**可写作左右值的组合：

$$
\mathbf U(x,t)=\sum_{i=I+1}^{m}{\alpha_i\mathbf K^{(i)}}+\sum_{i=1}^{I}{\beta_i\mathbf K^{(i)}},\quad \lambda_I<\frac{x}{t}<\lambda_{I+1}
$$

每个波的波强为：

$$
(\Delta \mathbf U)_i=(\beta_i - \alpha_i)\mathbf K^{(i)}
$$

![线性双曲系统 Riemann 问题的波系结构](576a0a96-ec40-4955-a4fc-161af6a132c0.png)
_线性双曲系统 Riemann 问题的波系结构_

此外还需理解**依赖域、决定域、影响域**的概念。

### 1.4 守恒律

守恒律的四种等价积分形式（控制域 $$[x_L,x_R]\times[t_1,t_2]$$）：

**形式 I**（对空间积分后求时间导数）：

$$
\frac{\text d}{\text dt} \int_{x_L}^{x_R} \mathbf{U}(x, t) \, dx = \mathbf{F}(\mathbf{U}(x_L, t)) - \mathbf{F}(\mathbf{U}(x_R, t))
$$

**形式 II**（再对时间积分）：

$$
\int_{x_L}^{x_R} \mathbf{U}(x, t_2)dx = \int_{x_L}^{x_R} \mathbf{U}(x, t_1)dx + \int_{t_1}^{t_2} \mathbf{F}(\mathbf{U}(x_L, t))dt - \int_{t_1}^{t_2} \mathbf{F}(\mathbf{U}(x_R, t)) \, dt
$$

**形式 III**（Green 定理转化为线积分）：

$$
\oint[\mathbf Udx-\mathbf {F(U)}dt]=\mathbf 0
$$

**形式 IV**（弱解形式，引入试函数 $$\phi$$，允许间断）：

$$
\int_{0}^{+\infty} \int_{-\infty}^{+\infty} \left[ \phi_t \mathbf{U} + \phi_x \mathbf{F}(\mathbf{U}) \right] \, dx \, dt = - \int_{-\infty}^{+\infty} \phi(x, 0) \mathbf{U}(x, 0) \, dx
$$

### 1.5 非线性波的演化

对无粘非线性对流方程 $$u_t+\lambda(u)u_x=0,\;\lambda(u)=f'(u)$$，沿特征线解 $$u$$ 不变（特征线仍为直线），但不同 $$u$$ 处特征线斜率不同，可能导致**间断的形成**。

![特征线相交导致间断形成](4e22e016-3ee0-46ac-ab14-5cf4195c8afb.png)
_特征线相交导致间断形成_

物理上需要 $$u_{xx}$$ 项对抗波变陡峭，或允许间断形成**激波**。

**Rankine-Hugoniot 条件**：

$$
\Delta f = S\Delta u
$$

**熵条件**：

$$
\lambda(u_L)>S>\lambda(u_R)
$$

不满足熵条件的激波（稀疏激波）结构不稳定，实际对应**稀疏波**。

> 经典示例——无粘 Burgers 方程的 Riemann 问题
{: .prompt-info }

对 $$u_t+(\frac{u^2}{2})_x=0$$：

- 当 $$u_L > u_R$$，形成激波：$$S=\frac{1}{2}(u_L+u_R)$$
- 当 $$u_L \le u_R$$，形成稀疏波：$$u(x,t)=x/t$$ 在 $$u_L < x/t < u_R$$ 区域

### 1.6 特征场分类与广义黎曼不变量

将标量方程的结果推广到双曲系统守恒律。

**特征场分类**：

- **线性退化**：若 $$\nabla \lambda_i(\mathbf U)\cdot\mathbf K^{(i)}(\mathbf U)=0,\;\forall \mathbf U\in \mathbb R^m$$
- **本质非线性**：若 $$\nabla \lambda_i(\mathbf U)\cdot\mathbf K^{(i)}(\mathbf U)\neq 0,\;\forall \mathbf U\in \mathbb R^m$$

**广义黎曼不变量**：考虑第 $$i$$ 族简单波区（所有状态变量通过单一标量参数 $$\xi(x,t)$$ 表达），满足：

$$
\frac{\text dw_1}{k_1^{(i)}}=\frac{\text dw_2}{k_2^{(i)}}=
\frac{\text dw_3}{k_3^{(i)}}=\dots=\frac{\text dw_m}{k_m^{(i)}}
$$

该关系表明跨越 $$\lambda_i$$ 简单波的两状态在状态空间中沿与 $$\mathbf K^{(i)}$$ 相切的积分曲线相连。沿该曲线存在 $$m-1$$ 个独立的**黎曼不变量**（保持常数），剩余一个参数沿曲线变化。

### 1.7 简单波解小结

![Riemann 问题的简单波解结构](d6e96101-b833-421b-bd64-655b0f1dc936.png)
_Riemann 问题的简单波解结构_

| 波类型 | 条件 | 关键关系 |
|--------|------|----------|
| 激波 | R-H 条件 + 熵条件 | $$\Delta\mathbf F=S_i\Delta\mathbf U$$ |
| 接触间断 | 线性退化场 | R-H 条件 + 广义黎曼不变量 + $$\lambda_i(\mathbf U_L)=\lambda_i(\mathbf U_R)=S_i$$ |
| 膨胀波 | 本质非线性场 | 广义黎曼不变量 + $$\lambda_i(\mathbf U_L)<\lambda_i(\mathbf U_R)$$ |

> **附：双曲系统的四种解法**——以线化气体方程为例，可采用特征线法、特征变量法、Riemann 不变量法、以及数值方法进行求解。
{: .prompt-tip }

---

## 2. Euler 方程的性质

### 2.1 一维 Euler 方程

**守恒形式**：

$$
\mathbf U_t+\mathbf {F(U)}_x=\mathbf 0
$$

其中：

$$
\mathbf U = \begin{bmatrix}\rho \\ \rho u \\ E\end{bmatrix},\quad
\mathbf F = \begin{bmatrix}\rho u \\ \rho u^2 + p \\ u(E + p)\end{bmatrix}
$$

**拟线性形式** $$\mathbf U_t+\mathbf A(\mathbf U)\mathbf U_x=0$$ 中：

$$
\mathbf{A(U)}=\begin{bmatrix}
0 & 1 & 0 \\
\frac{1}{2}(\gamma-3)u^2 & (3-\gamma)u & \gamma-1 \\
(\gamma-1)u^3-\gamma uE/\rho & \gamma E/\rho-\frac{3}{2}(\gamma-1)u^2 & \gamma u
\end{bmatrix}
$$

**特征值与右特征向量**：

$$
\lambda_1=u-a,\quad \lambda_2=u,\quad \lambda_3=u+a
$$

$$
\mathbf{K}^{(1)}=\begin{bmatrix}1 \\ u-a \\ H-ua\end{bmatrix},\quad
\mathbf{K}^{(2)}=\begin{bmatrix}1 \\ u \\ \frac{1}{2}u^2\end{bmatrix},\quad
\mathbf{K}^{(3)}=\begin{bmatrix}1 \\ u+a \\ H+ua\end{bmatrix}
$$

通量满足**齐次性**：$$\mathbf{F(U)=A(U)U}$$。

**原始变量形式**（$$\mathbf{W}=(\rho,u,p)^T$$）：
- 特征值不变
- 在光滑解下与守恒形式等价，但**在间断处给出错误解**
- 右特征向量方向不同（可数乘放缩），左特征向量可由 $$\mathbf L=\mathbf R^{-1}$$ 计算

沿特征线 $$\text dx/\text dt=\lambda_i$$ 的特征关系：

$$
\begin{align*}
&dp-\rho a\,du=0,\quad \text dx/\text dt=\lambda_1=u-a\\
&dp-a^2 d\rho=0,\quad \text dx/\text dt=\lambda_2=u\\
&dp+\rho a\,du=0,\quad \text dx/\text dt=\lambda_3=u+a
\end{align*}
$$

若用熵 $$s=c_v\ln(p/\rho^\gamma)+C_0$$ 替代 $$p$$，则在光滑区域熵沿流线不变：$$s_t+us_x=0$$。

### 2.2 一维 Riemann 问题的波系结构

三个特征值对应三个波，将空间分为四个常值状态：

![一维 Euler 方程 Riemann 问题解的结构](1724919424634.png)
_一维 Euler 方程 Riemann 问题解的结构_

- **接触间断**（$$\lambda_2=u$$）：由广义黎曼不变量得 $$p=\text{constant}$$，$$u=\text{constant}$$，仅密度跳跃
- **稀疏波**（$$\lambda_1=u-a,\;\lambda_3=u+a$$）：Riemann 不变量为 $$I_1=u\pm\frac{2a}{\gamma-1}$$，$$I_2=s$$
- **激波**（$$\lambda_1$$ 或 $$\lambda_3$$）：通过 R-H 条件连接两侧

> 激波关系详见空气动力学教材，核心是变换到激波静止参考系后应用 R-H 条件。
{: .prompt-info }

### 2.3 多维 Euler 方程

**二维守恒形式** $$\mathbf U_t + \mathbf F_x + \mathbf G_y = 0$$：

$$
\mathbf U =\begin{bmatrix}\rho \\ \rho u \\ \rho v \\ \rho E\end{bmatrix},\quad
\mathbf F =\begin{bmatrix}\rho u \\ \rho u^2+p \\ \rho uv \\ (\rho E+p)u\end{bmatrix},\quad
\mathbf G =\begin{bmatrix}\rho v \\ \rho uv \\ \rho v^2+p \\ (\rho E+p)v\end{bmatrix}
$$

其中 $$E=e+\frac{1}{2}(u^2+v^2)$$，$$H=E+p/\rho$$。

$$\mathbf A=\partial\mathbf F/\partial\mathbf U$$ 的特征值为 $$\lambda=\{u-a,u,u,u+a\}$$，对应特征场：
- 1, 4 波：本质非线性（激波或稀疏波）
- 2, 3 波：线性退化（密度间断与切向速度间断）

**旋转不变性**：对任意 $$\theta$$，

$$
\cos \theta\mathbf{F(U)}+\sin \theta\mathbf{G(U)}=\mathbf T^{-1}\mathbf{F(TU)}
$$

其中 $$\mathbf T(\theta)$$ 为旋转矩阵。由此可证 Euler 方程在时间方向上的双曲性。

![二维控制体积](1724921881287.png)
_二维控制体积_

**三维 Euler 方程**具有类似结构：5 个特征值 $$\{u-a,u,u,u,u+a\}$$，同样满足旋转不变性。

**x-方向分裂的 Riemann 问题**中，切向速度 $$v,w$$ 仅在中间波（剪切波）处跳跃：

![三维 x-方向分裂 Riemann 问题](QQ_1781187295702.png)
_三维 x-方向分裂 Riemann 问题的波系结构_

> 剪切波导致左右两侧切向速度跳跃是近似 Riemann 求解器的常见问题来源。
{: .prompt-warning }

### 2.4 守恒形式 vs 非守恒形式

虽然守恒和非守恒形式在光滑区等价，但**非守恒形式对激波等间断会给出错误解**。即使非守恒形式的方程在数学上也可写成某种"守恒形式"，这种守恒不具备物理意义。

以浅水波方程为例：守恒形式的右行激波波速为 $$S = u_R+Q/\phi_R$$，而原始变量形式给出 $$\hat S = u_R+\hat Q/\phi_R$$，一般有 $$\hat S\le S$$——非守恒形式低估了激波速度。

---

## 3. Euler 方程的 Riemann 问题精确解

### 3.1 求解策略

一维 Euler 方程 Riemann 问题：

$$
\mathbf U(x,0)=
\begin{cases}
\mathbf U_L  & \text{if } x<0 \\
\mathbf U_R  & \text{if } x>0
\end{cases}
$$

解的结构包含三个波，中间为**星区**（接触间断两侧分别为 $$\rho_{*L},\rho_{*R}$$，但共用 $$p_*,u_*$$）：

![Riemann 问题精确解结构](8ceaf704-81bf-46a4-9844-47166fb724ba.png)
_Riemann 问题精确解结构——主要求解量：$$p_*,u_*,\rho_{*L},\rho_{*R}$$_

### 3.2 压力与速度方程

星区压力 $$p_*$$ 是下述方程的解：

$$
f(p, \mathbf{W}_L, \mathbf{W}_R) \equiv f_L(p, \mathbf{W}_L) + f_R(p, \mathbf{W}_R) + \Delta u = 0,\quad \Delta u \equiv u_R - u_L
$$

其中 $$f_L,f_R$$ 为穿越波系的速度差：

$$
f_L(p, \mathbf{W}_L) = \begin{cases}
(p - p_L) \left[ \dfrac{A_L}{p + B_L} \right]^{\frac{1}{2}} & \text{if } p > p_L \text{ (激波)}, \\[10pt]
\dfrac{2a_L}{(\gamma - 1)} \left[ \left( \dfrac{p}{p_L} \right)^{\frac{\gamma - 1}{2\gamma}} - 1 \right] & \text{if } p \leq p_L \text{ (稀疏波)},
\end{cases}
$$

$$
f_R(p, \mathbf{W}_R) = \begin{cases}
(p - p_R) \left[ \dfrac{A_R}{p + B_R} \right]^{\frac{1}{2}} & \text{if } p > p_R \text{ (激波)}, \\[10pt]
\dfrac{2a_R}{(\gamma - 1)} \left[ \left( \dfrac{p}{p_R} \right)^{\frac{\gamma - 1}{2\gamma}} - 1 \right] & \text{if } p \leq p_R \text{ (稀疏波)},
\end{cases}
$$

常数：

$$
A_L = \frac{2}{(\gamma+1)\rho_L},\; B_L = \frac{(\gamma-1)}{(\gamma+1)}p_L,\;
A_R = \frac{2}{(\gamma+1)\rho_R},\; B_R = \frac{(\gamma-1)}{(\gamma+1)}p_R
$$

解得 $$p_*$$ 后，星区速度：

$$
u_*=\frac{1}{2}(u_L+u_R)+\frac{1}{2}[f_R(p_*)-f_L(p_*)]
$$

再由波系穿越关系求 $$\rho_{*L}$$ 和 $$\rho_{*R}$$。

> **推导要点**：激波分支通过参考系变换和 R-H 关系得到质量流量表达式，再代入激波密度比关系消去中间变量。稀疏波分支利用等熵关系与广义黎曼不变量直接求解。
{: .prompt-tip }

### 3.3 压力函数的数值求解

$$f(p)$$ 的性质：

- $$f'>0$$：单调递增，保证解唯一
- $$f''<0$$：凹函数

**压力正性条件**（防止真空）：

$$
(\Delta u)_\text{crit}\equiv \frac{2a_L}{\gamma-1}+\frac{2a_R}{\gamma-1}>u_R-u_L
$$

![不同波系配置下压力函数的行为](17d35143-7389-48aa-af49-1c47ee6785be.png)
_不同 $$\Delta u$$ 下压力函数的行为：$$(\Delta u)_3<(\Delta u)_2<(\Delta u)_1$$_

通过 $$f(p_{\min})$$ 和 $$f(p_{\max})$$ 可判断波系类型。当 $$\Delta u$$ 过大时可能得到负解，对应真空。

### 3.4 含真空的 Riemann 问题

真空指 $$\rho=0,p=0$$，但速度 $$u_0$$ 可为有限值。包含两种情况：

- 初始状态含真空（Case 1 / Case 2）
- 非真空状态发展出真空区（Case 3）

关键性质：
- 激波**不能**与真空区相邻
- 接触间断**可以**与真空区相邻
- 真空区速度 $$u_0$$ 即物质与真空边界的速度

**Case 1**（右真空）：稀疏波 + 接触间断，波尾处即为接触间断，波速 $$S_{*L}\equiv u_0=u_L+\frac{2a_L}{\gamma-1}$$。

![右真空情形](77554835-5d20-4aab-b535-bd27c7c0a4f0.png)
_右真空 Riemann 问题解结构_

**Case 3**（产生真空）：由 $$S_{*L}\le S_{*R}$$ 判定，等价于前述压力正性条件的违反。

![产生真空情形](5541297b-dd88-4b88-a940-19eda6ee089d.png)
_左右稀疏波之间产生真空区_

### 3.5 分裂多维 Riemann 问题

对三维 x-方向的 Riemann 问题，切向速度 $$v,w$$ 仅在中间波（接触间断/剪切波）处间断。对任意被动标量 $$q$$ 满足对流方程 $$q_t+uq_x+vq_y+wq_z=0$$，其在 Riemann 问题中的解为：

$$
q(x,t)=\begin{cases}
q_L & \text{if}\quad \dfrac{x}{t}<u_* \\[8pt]
q_R & \text{if}\quad \dfrac{x}{t}>u_*
\end{cases}
$$

---

## 总结

本文覆盖了该书前三章核心理论部分（第 2–4 章），建立了从双曲型 PDE 数学基础到 Euler 方程 Riemann 问题精确解的完整知识链。后续章节将进入数值方法部分，包括 Godunov 方法、近似 Riemann 求解器和高阶格式等内容。
