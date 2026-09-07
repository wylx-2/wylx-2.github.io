---
title: 《Riemann Solvers》读书笔记 — 第 3 章：Euler 方程的性质
author: wylx2
date: 2026-06-30 10:00:00 +0800
categories: [读书笔记, 《Riemann Solvers and Numerical Methods for Fluid Dynamics》]
tags: [读书笔记, CFD, Riemann, 流体力学, Euler方程]
description: 《Riemann Solvers and Numerical Methods for Fluid Dynamics》第 3 章笔记：Euler 方程的性质，涵盖一维与多维守恒形式、特征结构、简单波解以及守恒与非守恒形式的对比。
math: true
media_subpath: /assets/img/riemann-solvers/
---

> 本文是 E.F. Toro 所著 *Riemann Solvers and Numerical Methods for Fluid Dynamics* 读书笔记系列的一篇，对应原书第 3 章。
> 系列笔记：[第 2 章 · 双曲型偏微分方程基本概念](/posts/riemann-solvers-ch2-hyperbolic-pdes/) · [第 3 章 · Euler 方程的性质](/posts/riemann-solvers-ch3-euler-equations/) · [第 4 章 · Euler 方程的 Riemann 问题精确解](/posts/riemann-solvers-ch4-riemann-problem/)。原书第 1 章为引言，第 5 章起的 Riemann 求解器与高阶格式内容待后续补充。
{: .prompt-info }

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
