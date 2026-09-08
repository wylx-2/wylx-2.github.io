---
title: 《Riemann Solvers》读书笔记 — 第 4 章：Euler 方程的 Riemann 问题精确解
author: wylx2
date: 2026-06-30 12:00:00 +0800
categories: [读书笔记, 《Riemann Solvers and Numerical Methods for Fluid Dynamics》]
tags: [读书笔记, CFD, Riemann, 流体力学, Euler方程]
description: 《Riemann Solvers and Numerical Methods for Fluid Dynamics》第 4 章笔记：Euler 方程 Riemann 问题的精确解法，包括星区压力方程、左右行波函数、压力的数值求解及含真空情形。
math: true
media_subpath: /assets/img/riemann-solvers/
---

> 本文是 E.F. Toro 所著 *Riemann Solvers and Numerical Methods for Fluid Dynamics* 读书笔记系列的一篇，对应原书第 4 章。
> 系列笔记：[第 2 章 · 双曲型偏微分方程基本概念](/posts/riemann-solvers-ch2-hyperbolic-pdes/) · [第 3 章 · Euler 方程的性质](/posts/riemann-solvers-ch3-euler-equations/) · [第 4 章 · Euler 方程的 Riemann 问题精确解](/posts/riemann-solvers-ch4-riemann-problem/) · [第 10 章 · HLL 与 HLLC 黎曼求解器](/posts/riemann-solvers-ch10-hll-hllc-solvers/)。原书第 1 章为引言，其余章节待后续补充。
{: .prompt-info }

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
