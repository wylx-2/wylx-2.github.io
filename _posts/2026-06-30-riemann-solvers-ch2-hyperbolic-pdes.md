---
title: 《Riemann Solvers》读书笔记 — 第 2 章：双曲型偏微分方程基本概念
author: wylx2
date: 2026-06-30 08:00:00 +0800
categories: [读书笔记, 《Riemann Solvers and Numerical Methods for Fluid Dynamics》]
tags: [读书笔记, CFD, Riemann, 流体力学, 双曲系统]
description: 《Riemann Solvers and Numerical Methods for Fluid Dynamics》第 2 章笔记：双曲型偏微分方程的基本概念，包括拟线性系统、特征线法、守恒律、Rankine-Hugoniot 条件、熵条件与广义黎曼不变量。
math: true
media_subpath: /assets/img/riemann-solvers/
---

> 本文是 E.F. Toro 所著 *Riemann Solvers and Numerical Methods for Fluid Dynamics* 读书笔记系列的一篇，对应原书第 2 章。
> 系列笔记：[第 2 章 · 双曲型偏微分方程基本概念](/posts/riemann-solvers-ch2-hyperbolic-pdes/) · [第 3 章 · Euler 方程的性质](/posts/riemann-solvers-ch3-euler-equations/) · [第 4 章 · Euler 方程的 Riemann 问题精确解](/posts/riemann-solvers-ch4-riemann-problem/)。原书第 1 章为引言，第 5 章起的 Riemann 求解器与高阶格式内容待后续补充。
{: .prompt-info }

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
