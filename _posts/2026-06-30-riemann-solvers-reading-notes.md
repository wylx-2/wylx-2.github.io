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

本文是 E.F. Toro 所著 *Riemann Solvers and Numerical Methods for Fluid Dynamics* 的读书笔记，目前涵盖第 2–4 章的核心内容，其余章节待后续补充。

---

# 1 Introduction — 引言（待补充）

> 原书第 1 章为引言，介绍计算流体力学背景、有限体积法基本概念以及 Riemann 问题在 Godunov 方法中的核心地位，后续补充。
{: .prompt-info }

---

# 2 Notions on Hyperbolic Partial Differential Equations — 双曲型偏微分方程基本概念

本章介绍双曲型偏微分方程的基本性质，为后续 Euler 方程和 Riemann 问题的讨论奠定数学基础。

## 2.1 Quasi–Linear Equations: Basic Concepts — 拟线性系统基本概念

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

此外求解系统需要指定初边值条件。

**守恒律**形式：

$$
\mathbf U_t+\mathbf F(\mathbf U)_x=0
$$

**Jacobi 矩阵**：

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

## 2.2 The Linear Advection Equation — 线性对流方程

考虑线性对流方程：

$$
u_t+au_x=0,\quad u(x,0)=u_0(x),\quad -\infty <x<\infty,\;t>0
$$

**特征线法**：将 $$u$$ 视为 $$t$$ 的函数 $$u=u(x(t),t)$$，则

$$
\frac{\text d u}{\text d t}=\frac{\partial u}{\partial t}+\frac{\text d x}{\text d t}\frac{\partial u}{\partial x}
$$

若特征线满足 $$\text dx/\text dt=a$$，则 $$\text du/\text dt=0$$，即 $$u$$ 沿特征线不变。于是解为：

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

## 2.3 Linear Hyperbolic Systems — 线性双曲系统

考察：

$$
\mathbf U_t+\mathbf {AU}_x = \mathbf 0
$$

**对角化**：$$\mathbf A = \mathbf {K\Lambda K}^{-1}$$

**特征变量**：$$\mathbf {W} = \mathbf K^{-1}\mathbf U$$，解耦得到：

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

## 2.4 Conservation Laws — 守恒律

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

### 2.4.1 非线性波的演化

对无粘非线性对流方程 $$u_t+\lambda(u)u_x=0,\;\lambda(u)=f'(u)$$，沿特征线 $$\text dx/\text dt=\lambda(u)$$ 解 $$u$$ 不变（特征线仍为直线），但不同 $$u$$ 处特征线斜率不同，可能导致**间断的形成**。

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

> **经典示例**——无粘 Burgers 方程的 Riemann 问题
{: .prompt-info }

对 $$u_t+(\frac{u^2}{2})_x=0$$：

- 当 $$u_L > u_R$$，形成激波：

$$
S=\frac{1}{2}(u_L+u_R),\quad
u(x,t)=\begin{cases}
u_L  & \text{if } x-St<0 \\
u_R  & \text{if } x-St>0
\end{cases}
$$

- 当 $$u_L \le u_R$$，形成稀疏波：

$$
u(x,t)=\begin{cases}
u_L  & \text{if } x/t\le u_L \\
x/t  & \text{if } u_L<x/t<u_R \\
u_R  & \text{if } x/t\ge u_R
\end{cases}
$$

### 2.4.2 特征场分类与广义黎曼不变量

将标量方程的结果推广到双曲系统守恒律。

**特征场分类**：

- **线性退化**（linearly degenerate）：若 $$\nabla \lambda_i(\mathbf U)\cdot\mathbf K^{(i)}(\mathbf U)=0,\;\forall \mathbf U\in \mathbb R^m$$
- **本质非线性**（genuinely nonlinear）：若 $$\nabla \lambda_i(\mathbf U)\cdot\mathbf K^{(i)}(\mathbf U)\neq 0,\;\forall \mathbf U\in \mathbb R^m$$

Rankine-Hugoniot 条件：

$$
\Delta \mathbf F = S_i\Delta \mathbf U
$$

其中间断波速 $$S_i$$ 与 $$\lambda_i$$ 特征场相关。对线性系统，$$S_i=\lambda_i$$。

**广义黎曼不变量**（generalized Riemann invariants）：

考虑拟线性双曲系统 $$\mathbf W_t + \mathbf A(\mathbf W) \mathbf W_x = 0$$。假定解处于一个只包含第 $$i$$ 族波的**简单波区**，所有状态变量通过一个标量参数 $$\xi = \xi(x,t)$$ 表达。第 $$i$$ 族波满足特征方程 $$\xi_t + \lambda_i \xi_x = 0$$，代入原方程化简得：

$$
\left[ \mathbf A(\mathbf W) - \lambda_i \mathbf I \right] \frac{d\mathbf W}{d\xi} = 0
$$

对应于右特征向量，即 $$\text d \mathbf W \parallel \mathbf K^{(i)}$$，分量形式为：

$$
\frac{\text dw_1}{k_1^{(i)}}=\frac{\text dw_2}{k_2^{(i)}}=
\frac{\text dw_3}{k_3^{(i)}}=\dots=\frac{\text dw_m}{k_m^{(i)}}
$$

该关系表明，跨越 $$\lambda_i$$ 简单波的两个状态在状态空间中沿着一条**积分曲线**相连，该曲线处处与右特征向量场 $$\mathbf K^{(i)}$$ 相切。沿此曲线，存在 $$m-1$$ 个独立的函数（称为**黎曼不变量**）保持为常数，剩余一个参数沿曲线变化。

### 2.4.3 简单波解小结

![Riemann 问题的简单波解结构](d6e96101-b833-421b-bd64-655b0f1dc936.png)
_Riemann 问题的简单波解结构_

| 波类型 | 条件 | 关键关系 |
|--------|------|----------|
| 激波 (shock) | R-H 条件 + 熵条件 | $$\Delta\mathbf F=S_i\Delta\mathbf U$$ |
| 接触间断 (contact) | 线性退化场 | R-H 条件 + 广义黎曼不变量 + $$\lambda_i(\mathbf U_L)=\lambda_i(\mathbf U_R)=S_i$$ |
| 膨胀波 (rarefaction) | 本质非线性场 | 广义黎曼不变量 + $$\lambda_i(\mathbf U_L)<\lambda_i(\mathbf U_R)$$ |

> **附：双曲系统的四种解法**——以线化气体方程为例，可采用特征线法、特征变量法、Riemann 不变量法以及数值方法进行求解。
{: .prompt-tip }

---

# 3 Some Properties of the Euler Equations — Euler 方程的性质

## 3.1 The One–Dimensional Euler Equations — 一维 Euler 方程

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

**非守恒形式**（原始变量 $$\mathbf{W}=(\rho,u,p)^T$$）：

$$
\mathbf{A(W)}=\begin{bmatrix}
u & \rho & 0 \\
0 & u & 1/\rho \\
0 & \rho a^2 & u
\end{bmatrix}
$$

特征值与守恒形式一致。在光滑解下非守恒形式与守恒形式等价，但**在间断处给出错误解**。

右特征向量：

$$
\mathbf{K}^{(1)}=\alpha_1\begin{bmatrix}1 \\ -a/\rho \\ a^2\end{bmatrix},\quad
\mathbf{K}^{(2)}=\alpha_2\begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix},\quad
\mathbf{K}^{(3)}=\alpha_3\begin{bmatrix}1 \\ a/\rho \\ a^2\end{bmatrix}
$$

沿特征线 $$\text dx/\text dt=\lambda_i$$ 的特征关系（$$\mathbf{L}^{(i)}\cdot d\mathbf W=0$$）：

$$
\begin{align*}
&dp-\rho a\,du=0,\quad \text dx/\text dt=\lambda_1=u-a\\
&dp-a^2 d\rho=0,\quad \text dx/\text dt=\lambda_2=u\\
&dp+\rho a\,du=0,\quad \text dx/\text dt=\lambda_3=u+a
\end{align*}
$$

若用熵 $$s=c_v\ln(p/\rho^\gamma)+C_0$$ 替代 $$p$$，则 $$\mathbf{W}=(\rho,u,s)^T$$ 下沿流线 $$s_t+us_x=0$$（光滑区域熵不变）。

### 3.1.1 Riemann 问题的简单波

三个特征值形成的三个波将空间分为四个常值状态。其中 $$\mathbf{K}^{(2)}$$ 特征空间关联的波是接触间断，而 $$\mathbf{K}^{(1)}$$ 和 $$\mathbf{K}^{(3)}$$ 为稀疏波或激波。

![一维 Euler 方程 Riemann 问题解的结构](1724919424634.png)
_在 x-t 平面上 Riemann 问题解的结构_

**接触间断**（contact discontinuity）：由广义黎曼不变量得两侧 $$p=\text{constant}$$，$$u=\text{constant}$$，仅密度跳跃。

**稀疏波**（rarefaction wave）：广义黎曼不变量：
- 对 $$\lambda_1=u-a$$：$$I_1=u+\dfrac{2a}{\gamma-1},\; I_2=s$$
- 对 $$\lambda_3=u+a$$：$$I_1=u-\dfrac{2a}{\gamma-1},\; I_2=s$$

波内 $$\rho,u,p$$ 连续光滑变化，具有扇形结构。

**激波**（shock wave）：通过 R-H 条件连接两侧，变换到激波静止参考系后求解。以 3-波系为例：

$$
\frac{\rho_*}{\rho_R} = \frac{(\gamma + 1)(M_R - M_S)^2}{(\gamma - 1)(M_R - M_S)^2 + 2},\quad
\frac{p_*}{p_R} = \frac{2\gamma(M_R - M_S)^2 - (\gamma - 1)}{(\gamma + 1)}
$$

其中 $$M_R=u_R/a_R,\; M_S=S_3/a_R$$。1-波系类似。

## 3.2 Multi–Dimensional Euler Equations — 多维 Euler 方程

### 3.2.1 Two–Dimensional Equations in Conservative Form — 二维守恒形式

$$
\mathbf U_t + \mathbf F_x + \mathbf G_y = 0
$$

其中：

$$
\mathbf U =\begin{bmatrix}\rho \\ \rho u \\ \rho v \\ \rho E\end{bmatrix},\quad
\mathbf F =\begin{bmatrix}\rho u \\ \rho u^2+p \\ \rho uv \\ (\rho E+p)u\end{bmatrix},\quad
\mathbf G =\begin{bmatrix}\rho v \\ \rho uv \\ \rho v^2+p \\ (\rho E+p)v\end{bmatrix}
$$

其中 $$E=e+\frac{1}{2}(u^2+v^2)$$，$$H=E+p/\rho$$。

$$\mathbf A=\partial\mathbf F/\partial\mathbf U$$ 的特征值为 $$\lambda=\{u-a,u,u,u+a\}$$，对应特征场：
- 1, 4 波：本质非线性（激波或稀疏波）
- 2, 3 波：线性退化（密度间断与切向速度间断）

**旋转不变性**：

![二维控制体积](1724921881287.png)

对任意 $$\theta$$，$$\cos \theta\mathbf{F(U)}+\sin \theta\mathbf{G(U)}=\mathbf T^{-1}\mathbf{F(TU)}$$，其中 $$\mathbf T(\theta)$$ 为旋转矩阵。

### 3.2.2 Three–Dimensional Equations in Conservative Form — 三维守恒形式

三维 Euler 方程具有 5 个特征值 $$\{u-a,u,u,u,u+a\}$$，同样满足旋转不变性。

**x-方向分裂 Riemann 问题**：

$$
\mathbf U_t+\mathbf {F(U)}_x=\mathbf 0,\quad
\mathbf U(x,0)=\begin{cases}
\mathbf U_L  & \text{if } x<0 \\
\mathbf U_R  & \text{if } x>0
\end{cases}
$$

其中 $$\mathbf U=[\rho,\rho u,\rho v,\rho w,E]^T$$。中间波包含剪切波，切向速度跳跃：

![三维 x-方向分裂 Riemann 问题](QQ_1781187295702.png)
_三维 x-方向分裂 Riemann 问题的波系结构_

> 剪切波导致左右两侧切向速度跳跃，是近似 Riemann 求解器的常见问题来源。
{: .prompt-warning }

## 3.3 Conservative Versus Non–Conservative Formulations — 守恒与非守恒形式

虽然守恒和非守恒形式在光滑区等价，但**非守恒形式对激波等间断会给出错误解**。即使非守恒形式的方程在数学上也可写成某种"守恒形式"，这种守恒不具备物理意义。

以浅水波方程为例，守恒形式的右行激波波速为 $$S = u_R+Q/\phi_R$$，而原始变量形式给出 $$\hat S = u_R+\hat Q/\phi_R$$，一般有 $$\hat S\le S$$——非守恒形式低估了激波速度。

---

# 4 The Riemann Problem for the Euler Equations — Euler 方程的 Riemann 问题精确解

本章介绍理想气体 Euler 方程的 Riemann 问题精确求解过程。

## 4.1 Solution Strategy — 求解策略

一维 Euler 方程 Riemann 问题：

$$
\mathbf U(x,0)=\begin{cases}
\mathbf U_L  & \text{if } x<0 \\
\mathbf U_R  & \text{if } x>0
\end{cases}
$$

解的结构包含三个波，中间为**星区**（star region）——接触间断两侧分别为 $$\rho_{*L},\rho_{*R}$$，但共用 $$p_*,u_*$$：

![Riemann 问题精确解结构](8ceaf704-81bf-46a4-9844-47166fb724ba.png)
_Riemann 问题精确解结构——主要求解量：$$p_*,u_*,\rho_{*L},\rho_{*R}$$_

## 4.2 Equations for Pressure and Particle Velocity — 压力与速度方程

星区压力 $$p_*$$ 是下述方程的解：

$$
f(p, \mathbf{W}_L, \mathbf{W}_R) \equiv f_L(p, \mathbf{W}_L) + f_R(p, \mathbf{W}_R) + \Delta u = 0,\quad \Delta u \equiv u_R - u_L
$$

其中 $$f_L,f_R$$ 为穿越波系的速度差：

$$
f_L(p, \mathbf{W}_L) = \begin{cases}
(p - p_L) \left[ \dfrac{A_L}{p + B_L} \right]^{\frac{1}{2}} & \text{if } p > p_L \text{ (shock)}, \\[10pt]
\dfrac{2a_L}{(\gamma - 1)} \left[ \left( \dfrac{p}{p_L} \right)^{\frac{\gamma - 1}{2\gamma}} - 1 \right] & \text{if } p \leq p_L \text{ (rarefaction)},
\end{cases}
$$

$$
f_R(p, \mathbf{W}_R) = \begin{cases}
(p - p_R) \left[ \dfrac{A_R}{p + B_R} \right]^{\frac{1}{2}} & \text{if } p > p_R \text{ (shock)}, \\[10pt]
\dfrac{2a_R}{(\gamma - 1)} \left[ \left( \dfrac{p}{p_R} \right)^{\frac{\gamma - 1}{2\gamma}} - 1 \right] & \text{if } p \leq p_R \text{ (rarefaction)},
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

### 4.2.1 Function $$f_L$$ for a Left Shock — 左行激波的 f_L 函数

通过参考系变换和 R-H 关系，定义质量流量 $$Q_L \equiv \rho_L \hat u_L = \rho_{*L}\hat u_*$$，利用激波密度比关系消去中间变量，最终得到：

$$
u_* = u_L - (p_*-p_L)\left[\frac{A_L}{p_*+B_L}\right]^{\frac{1}{2}}
$$

### 4.2.2 Function $$f_L$$ for Left Rarefaction — 左行稀疏波的 f_L 函数

由等熵关系和广义黎曼不变量：

$$
\rho_{*L} = \rho_{L} \left( \frac{p_{*}}{p_{L}} \right)^\frac{1}{\gamma},\quad
a_{*L} = a_{L} \left( \frac{p_{*}}{p_{L}} \right)^\frac{\gamma-1}{2\gamma}
$$

代入 $$u_L + \dfrac{2a_L}{\gamma-1} = u_* + \dfrac{2a_*}{\gamma-1}$$ 得：

$$
u_* = u_L - \frac{2a_L}{\gamma-1}\left[\left(\frac{p_*}{p_L}\right)^{\frac{\gamma-1}{2\gamma}}-1\right]
$$

右行波系（3-波）类似。

## 4.3 Numerical Solution for Pressure — 压力的数值求解

### 4.3.1 Behavior of the Pressure Function — 压力函数的性质

$$f(p)$$：$$f'>0$$（单调递增），$$f''<0$$（凹函数）。

**压力正性条件**（防止真空）：

$$
(\Delta u)_\text{crit} \equiv \frac{2a_L}{\gamma-1}+\frac{2a_R}{\gamma-1} > u_R-u_L
$$

![不同波系配置下压力函数的行为](17d35143-7389-48aa-af49-1c47ee6785be.png)
_压力函数的行为：$$(\Delta u)_3<(\Delta u)_2<(\Delta u)_1$$_

通过 $$f(p_{\min})$$ 和 $$f(p_{\max})$$ 可判断波系类型。当 $$\Delta u$$ 过大时可能得到负解，对应真空。

## 4.6 The Riemann Problem in the Presence of Vacuum — 含真空的 Riemann 问题

真空指 $$\rho=0,p=0$$，但速度 $$u_0$$ 可为有限值。包含：
- 初始状态含真空
- 非真空状态发展出真空区

关键性质：激波**不能**与真空区相邻；接触间断**可以**与真空区相邻。

### 4.6.1 Case 1: Vacuum Right State — 右真空

稀疏波 + 接触间断，波速 $$S_{*L}\equiv u_0 = u_L+\dfrac{2a_L}{\gamma-1}$$。

![右真空情形](77554835-5d20-4aab-b535-bd27c7c0a4f0.png)
_右真空 Riemann 问题解结构_

### 4.6.3 Case 3: Generation of Vacuum — 真空产生

条件：$$S_{*L}\le S_{*R}\;\Rightarrow\; \dfrac{2a_L}{\gamma-1}+\dfrac{2a_R}{\gamma-1}\le u_R-u_L$$

![产生真空情形](5541297b-dd88-4b88-a940-19eda6ee089d.png)
_左右稀疏波之间产生真空区_

## 4.8 The Split Multi–Dimensional Case — 分裂多维情形

切向速度 $$v,w$$ 仅在中间波处间断。对被动标量 $$q$$ 满足 $$q_t+uq_x+vq_y+wq_z=0$$：

$$
q(x,t)=\begin{cases}
q_L & \text{if}\quad \dfrac{x}{t}<u_* \\[8pt]
q_R & \text{if}\quad \dfrac{x}{t}>u_*
\end{cases}
$$

---

# 5 Riemann Solvers（待补充）

> 原书第 5 章起介绍各类 Riemann 求解器，包括 Godunov 方法、近似 Riemann 求解器（HLL、HLLC、Roe、Osher 等），待后续补充。
{: .prompt-info }

# 6 Further Chapters — 高阶格式与后续章节（待补充）

> 原书后续章节涵盖 MUSCL 格式、TVD 方法、非结构化网格等内容，待后续补充。
{: .prompt-info }