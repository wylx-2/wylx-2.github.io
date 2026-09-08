---
title: 《Riemann Solvers》读书笔记 — 第 10 章：HLL 与 HLLC 黎曼求解器
author: wylx2
date: 2026-09-08 12:00:00 +0800
categories: [读书笔记, 《Riemann Solvers and Numerical Methods for Fluid Dynamics》]
tags: [读书笔记, CFD, Riemann, 流体力学, Euler方程, HLL, HLLC]
description: 《Riemann Solvers and Numerical Methods for Fluid Dynamics》第 10 章笔记：HLL 与 HLLC 近似黎曼求解器——双波与三波模型的推导、中间状态与通量表达式、波速估计方法及 HLLC 通量的计算步骤。
math: true
media_subpath: /assets/img/riemann-solvers/
---

> 本文是 E.F. Toro 所著 *Riemann Solvers and Numerical Methods for Fluid Dynamics* 读书笔记系列的一篇，对应原书第 10 章；本章及之后的笔记采用节选翻译＋脚注方式。
> 系列笔记：[第 2 章 · 双曲型偏微分方程基本概念](/posts/riemann-solvers-ch2-hyperbolic-pdes/) · [第 3 章 · Euler 方程的性质](/posts/riemann-solvers-ch3-euler-equations/) · [第 4 章 · Euler 方程的 Riemann 问题精确解](/posts/riemann-solvers-ch4-riemann-problem/) · [第 10 章 · HLL 与 HLLC 黎曼求解器](/posts/riemann-solvers-ch10-hll-hllc-solvers/)。原书第 1 章为引言，其余章节待后续补充。
{: .prompt-info }

Harten、Lax 与 van Leer（HLL）于 1983 年提出的近似黎曼求解器需要对由界面处初始间断发出的最快信号速度作出估计，从而形成一种以两道波刻画精确解结构的双波模型。更为精确的方法是 Toro 及其合作者于 1992 年引入的 HLLC：该方法采用三波模型，从而对中间波具有更高的分辨能力。

## 10.1 引言

为计算 Godunov 通量，Harten、Lax 与 van Leer 提出了一种近似求解黎曼问题的新途径，其核心思想是**为解假设一种由两道波分隔三个常数状态的波系构形**。
Davis 与 Einfeldt 二人独立地提出了多种计算波速的方法，使这一方法成为实用的格式，被称为 HLLE 黎曼求解器，或者 HLL 求解器。

然而，此类格式存在一个问题，其假设仅对含两个方程的双曲系统成立，例如一维浅水方程。对更大的系统，诸如 Euler 方程或分裂后的二维浅水方程，两波假设并不正确，会导致接触面、剪切波与物质界面等物理特征可能不能被精确分辨；当这些特征相对网格驻定时，由此产生的数值误差是不可接受的。为了解决这一问题，Einfeldt 对 HLLE 格式提出了一个称为 HLLM 的修正：采用一个线性分布代替 HLL 中的单一中间状态。此外，Toro、Spruce 与 Speares 于 1992 年提出了另一种方案：HLLC 黎曼求解器。HLLC 是一种三波模型，在黎曼问题解扇区的中间区域给出两个星区状态。

由于 Euler 方程在一维到三维空间中均只有三个不同的特征场，因此 HLLC 对 Euler 方程而言是一个完整的黎曼求解器——即 HLLC 的近似波结构包含精确问题的全部特征场。然而，对特征结构包含三个以上不同特征场的系统，HLLC 便不再完整，黎曼求解器的不完整性会损害中间波的分辨率，尤其是当这些波相对网格缓慢移动之时。因此，改进 HLLC 的一个显然途径，就是为所研究的系统接纳正确数目的特征场。这方面工作和其他发展有很多，这里不再列举。

## 10.2 黎曼问题（The Riemann Problem）

![图 10.1 x–t 平面上的控制体](fig10.1.png){: w="700" }
_**图 10.1**　x–t 平面上的控制体 $[x_L, x_R] \times [0, T]$。$S_L$ 与 $S_R$ 为黎曼问题之解所产生的最快信号速度_

略去黎曼问题的说明，考虑黎曼问题精确解所产生的全部波结构都包含于控制体 $[x_L, x_R] \times [0, T]$ 之内，如图 10.1 所示，即

$$
x_L \le TS_L,\qquad x_R \ge TS_R,
\tag{10.1}
$$

其中 $S_L$ 与 $S_R$ 分别为扰动初始数据状态 $\mathbf{U}_L$ 与 $\mathbf{U}_R$ 的最快信号速度，$T$ 为某一选定的时刻，于是给出积分形式

$$
\int_{x_L}^{x_R} \mathbf{U}(x,T)\,dx = \int_{x_L}^{x_R} \mathbf{U}(x,0)\,dx + \int_0^T \mathbf{F}(\mathbf{U}(x_L,t))\,dt - \int_0^T \mathbf{F}(\mathbf{U}(x_R,t))\,dt .
\tag{10.2}
$$

守恒律积分形式的细节见[本系列第 2 章](/posts/riemann-solvers-ch2-hyperbolic-pdes/)。计算上式右端，得

$$
\int_{x_L}^{x_R} \mathbf{U}(x,T)\,dx = x_R\mathbf{U}_R - x_L\mathbf{U}_L + T(\mathbf{F}_L - \mathbf{F}_R),
\tag{10.3}
$$

其中 $\mathbf{F}_L = \mathbf{F}(\mathbf{U}_L)$，$\mathbf{F}_R = \mathbf{F}(\mathbf{U}_R)$。我们称积分关系 (10.3) 为**相容性条件**。现在把 (10.2) 左端的积分拆为三个积分：

$$
\int_{x_L}^{x_R} \mathbf{U}(x,T)\,dx = \int_{x_L}^{TS_L} \mathbf{U}(x,T)\,dx + \int_{TS_L}^{TS_R} \mathbf{U}(x,T)\,dx + \int_{TS_R}^{x_R} \mathbf{U}(x,T)\,dx ,
$$

并计算右端第一、第三项，得

$$
\int_{x_L}^{x_R} \mathbf{U}(x,T)\,dx = \int_{TS_L}^{TS_R} \mathbf{U}(x,T)\,dx + (TS_L - x_L)\mathbf{U}_L + (x_R - TS_R)\mathbf{U}_R .
\tag{10.4}
$$

将 (10.4) 与 (10.3) 相比较，得

$$
\int_{TS_L}^{TS_R} \mathbf{U}(x,T)\,dx = T\left(S_R\mathbf{U}_R - S_L\mathbf{U}_L + \mathbf{F}_L - \mathbf{F}_R\right).
\tag{10.5}
$$

两边同时除以长度 $T(S_R - S_L)$——即 $T$ 时刻黎曼问题之解的波系在最慢与最快信号之间的宽度，得

$$
\frac{1}{T(S_R - S_L)}\int_{TS_L}^{TS_R} \mathbf{U}(x,T)\,dx = \frac{S_R\mathbf{U}_R - S_L\mathbf{U}_L + \mathbf{F}_L - \mathbf{F}_R}{S_R - S_L}.
\tag{10.6}
$$

于是，只要信号速度 $S_L$ 与 $S_R$ 已知，$T$ 时刻最慢与最快信号之间黎曼问题精确解的积分平均便是一个**已知的常数**；该常数即 (10.6) 的右端，记作

$$
\mathbf{U}_{hll} = \frac{S_R\mathbf{U}_R - S_L\mathbf{U}_L + \mathbf{F}_L - \mathbf{F}_R}{S_R - S_L}.
\tag{10.7}
$$

现在把守恒律的积分形式应用于图 10.1 的左半部分，即控制体 $[x_L, 0] \times [0, T]$，得

$$
\int_{TS_L}^{0} \mathbf{U}(x,T)\,dx = -TS_L\mathbf{U}_L + T(\mathbf{F}_L - \mathbf{F}_{0L}),
\tag{10.8}
$$

其中 $\mathbf{F}_{0L}$ 为沿 $t$ 轴的通量 $\mathbf{F}(\mathbf{U})$。解出 $\mathbf{F}_{0L}$，得

$$
\mathbf{F}_{0L} = \mathbf{F}_L - S_L\mathbf{U}_L - \frac{1}{T}\int_{TS_L}^{0} \mathbf{U}(x,T)\,dx .
\tag{10.9}
$$

在控制体 $[0, x_R] \times [0, T]$ 上求守恒律积分形式之值，则得

$$
\mathbf{F}_{0R} = \mathbf{F}_R - S_R\mathbf{U}_R + \frac{1}{T}\int_{0}^{TS_R} \mathbf{U}(x,T)\,dx .
\tag{10.10}
$$

读者容易验证，等式 $\mathbf{F}_{0L} = \mathbf{F}_{0R}$ 恰好给出相容性条件 (10.3)。由于假设的是黎曼问题的精确解，迄今所得的一切关系都是**精确的**。

## 10.3 HLL 近似黎曼求解器（The HLL Approximate Riemann Solver）

Harten、Lax 与 van Leer 提出了如下的近似黎曼求解器

$$
\tilde{\mathbf{U}}(x,t) = \begin{cases}
\mathbf{U}_L & \text{if } x/t \le S_L,\\
\mathbf{U}_{hll} & \text{if } S_L \le x/t \le S_R,\\
\mathbf{U}_R & \text{if } x/t \ge S_R,
\end{cases}
\tag{10.11}
$$

其中 $\mathbf{U}_{hll}$ 为由 (10.7) 给出的常数状态向量，速度 $S_L$ 与 $S_R$ 视为已知。图 10.2 给出了这一近似解的结构，称为 HLL 黎曼求解器。注意，该近似仅由两道波分隔的**三个常数状态**构成：星区只含单一常数状态，所有由中间波分隔的中间状态都被并拢为单一状态 $\mathbf{U}_{hll}$。相应的沿 $t$ 轴的通量 $\mathbf{F}_{hll}$ 由关系 (10.9) 或 (10.10) 求得，其中精确的被积函数换为近似解 (10.11)。注意，我们并**不**取 $\mathbf{F}_{hll} = \mathbf{F}(\mathbf{U}_{hll})$。

![图 10.2 HLL 近似黎曼求解器的波系结构](fig10.2.png){: w="700" }
_**图 10.2**　近似 HLL 黎曼求解器。星区内的解由单一状态 $\mathbf{U}_{hll}$ 构成，它与两个数据状态之间由速度为 $S_L$ 与 $S_R$ 的两道波分隔。_

有实际意义的非平凡情形是亚声速情形 $S_L \le 0 \le S_R$。把 (10.9) 或 (10.10) 中的被积函数代之以 (10.7) 的 $\mathbf{U}_{hll}$，得

$$
\mathbf{F}_{hll} = \mathbf{F}_L + S_L(\mathbf{U}_{hll} - \mathbf{U}_L),
$$

或

$$
\mathbf{F}_{hll} = \mathbf{F}_R + S_R(\mathbf{U}_{hll} - \mathbf{U}_R).
$$

把 (10.7) 代入上式，得 HLL 通量

$$
\mathbf{F}_{hll} = \frac{S_R\mathbf{F}_L - S_L\mathbf{F}_R + S_LS_R(\mathbf{U}_R - \mathbf{U}_L)}{S_R - S_L}.
\tag{10.12}
$$

于是，近似 Godunov 方法所用的 HLL 界面通量

$$
\mathbf{F}_{i+\frac{1}{2}}^{hll} = \begin{cases}
\mathbf{F}_L & \text{if } 0 \le S_L,\\[0.5ex]
\dfrac{S_R\mathbf{F}_L - S_L\mathbf{F}_R + S_LS_R(\mathbf{U}_R - \mathbf{U}_L)}{S_R - S_L} & \text{if } S_L \le 0 \le S_R,\\[1.5ex]
\mathbf{F}_R & \text{if } 0 \ge S_R .
\end{cases}
\tag{10.13}
$$

给定一个计算速度 $S_L$ 与 $S_R$ 的算法，就能通过 (10.13) 计算界面通量，从而构成一个近似 Godunov 方法。估计波速 $S_L$ 与 $S_R$ 的若干方法见第 10.5 节。

Harten、Lax 与 van Leer 证明了：Godunov 格式若收敛，则收敛于守恒律的弱解。事实上，他们证明的是，该近似解与守恒律的积分形式相容，也就是说若把精确解 $\mathbf{U}(x,t)$ 代入相容性条件 (10.3)，右端是保持不变的。

HLL 格式的一个缺点，在积分 (10.6) 中起作用的只是横跨波结构的**平均**，而完全没有考虑星区内黎曼问题解的空间变化。这一缺陷可以通过补回缺失的波而得到矫正。于是，Toro、Spruce 与 Speares 提出了所谓 HLLC 格式。

## 10.4 HLLC 近似黎曼求解器（The HLLC Approximate Riemann Solver）

### 10.4.1 若干有用关系

![图 10.3 HLLC 近似黎曼求解器的波系结构](fig10.3.png){: w="700" }
_**图 10.3**　HLLC 近似黎曼求解器。星区内的解由两个常数状态构成，二者之间由一道速度为 $S_*$ 的中间波分隔。_

补充一些积分关系，在控制体中再纳入一道速度为 $S_*$ 的中间波，如图 10.3 所示。把积分 (10.6) 的左端拆为两项，得

$$
\frac{1}{T(S_R - S_L)}\int_{TS_L}^{TS_R} \mathbf{U}(x,T)\,dx = \frac{1}{T(S_R - S_L)}\int_{TS_L}^{TS_*} \mathbf{U}(x,T)\,dx + \frac{1}{T(S_R - S_L)}\int_{TS_*}^{TS_R} \mathbf{U}(x,T)\,dx .
\tag{10.14}
$$

定义积分平均

$$
\mathbf{U}_L^* = \frac{1}{T(S_* - S_L)}\int_{TS_L}^{TS_*} \mathbf{U}(x,T)\,dx,\\
\mathbf{U}_R^* = \frac{1}{T(S_R - S_*)}\int_{TS_*}^{TS_R} \mathbf{U}(x,T)\,dx .
\tag{10.15}
$$

把 (10.15) 代入 (10.14) 并利用 (10.6)，相容性条件 (10.3) 成为

$$
\frac{S_* - S_L}{S_R - S_L}\mathbf{U}_L^* + \frac{S_R - S_*}{S_R - S_L}\mathbf{U}_R^* = \mathbf{U}_{hll},
\tag{10.16}
$$

其中 $\mathbf{U}_{hll}$ 由 (10.7) 给出。HLLC 近似黎曼求解器给出如下：

$$
\tilde{\mathbf{U}}(x,t) = \begin{cases}
\mathbf{U}_L & \text{if } x/t \le S_L,\\
\mathbf{U}_L^* & \text{if } S_L \le x/t \le S_*,\\
\mathbf{U}_R^* & \text{if } S_* \le x/t \le S_R,\\
\mathbf{U}_R & \text{if } x/t \ge S_R .
\end{cases}
\tag{10.17}
$$

我们寻求与之相应的 HLLC 数值通量，其定义为

$$
\mathbf{F}_{i+\frac{1}{2}}^{hllc} = \begin{cases}
\mathbf{F}_L & \text{if } 0 \le S_L,\\
\mathbf{F}_L^* & \text{if } S_L \le 0 \le S_*,\\
\mathbf{F}_R^* & \text{if } S_* \le 0 \le S_R,\\
\mathbf{F}_R & \text{if } 0 \ge S_R ,
\end{cases}
\tag{10.18}
$$

其中中间通量 $\mathbf{F}_L^*$ 与 $\mathbf{F}_R^*$ 尚待确定。图 10.3 给出了 HLLC 近似黎曼求解器的结构。对波速为 $S_L$、$S_*$、$S_R$ 的每一道波应用 Rankine–Hugoniot 关系，可得

$$
\mathbf{F}_L^* = \mathbf{F}_L + S_L(\mathbf{U}_L^* - \mathbf{U}_L),
\tag{10.19}
$$

$$
\mathbf{F}_R^* = \mathbf{F}_L^* + S_*(\mathbf{U}_R^* - \mathbf{U}_L^*),
\tag{10.20}
$$

$$
\mathbf{F}_R^* = \mathbf{F}_R + S_R(\mathbf{U}_R^* - \mathbf{U}_R).
\tag{10.21}
$$

三个方程，含四个未知向量 $\mathbf{U}_L^*$、$\mathbf{F}_L^*$、$\mathbf{U}_R^*$、$\mathbf{F}_R^*$。

### 10.4.2 Euler 方程的 HLLC 通量（The HLLC Flux for the Euler Equations）

为了求两个未知的中间通量 $\mathbf{F}_L^*$ 与 $\mathbf{F}_R^*$，只需求两个中间状态向量 $\mathbf{U}_L^*$ 与 $\mathbf{U}_R^*$。方程数少于未知向量数，需要附加条件。

关于中间波，精确解具有如下性质。
$$
p_L^* = p_R^* = p^*,\qquad u_L^* = u_R^* = u^*, \\
v_L^* = v_L,\quad v_R^* = v_R,\quad w_L^* = w_L,\quad w_R^* = w_R .\\
S_* = u^*
\tag{10.22}
$$

于是主要估计 $S_*$。现在，方程 (10.19) 与 (10.21) 可以改写为

$$
S_L\mathbf{U}_L^* - \mathbf{F}_L^* = S_L\mathbf{U}_L - \mathbf{F}_L,
\tag{10.23}
$$

$$
S_R\mathbf{U}_R^* - \mathbf{F}_R^* = S_R\mathbf{U}_R - \mathbf{F}_R,
\tag{10.24}
$$

其中 (10.23) 与 (10.24) 的右端均为已知的常数向量。我们还注意到 $\mathbf{U}$ 与 $\mathbf{F}$ 之间的有用关系，即

$$
\mathbf{F}(\mathbf{U}) = u\mathbf{U} + p\mathbf{D},\qquad \mathbf{D} = [0, 1, 0, 0, u]^T .
\tag{10.25}
$$

若波速 $S_L$ 与 $S_R$ 已知，对 (10.23)–(10.24) 的第一与第二个分量作代数运算，可得两个星区内压力的如下解：

$$
p_L^* = p_L + \rho_L(S_L - u_L)(S_* - u_L),\\
p_R^* = p_R + \rho_R(S_R - u_R)(S_* - u_R).
\tag{10.26}
$$

由 (10.22) 有 $p_L^* = p_R^*$；于是从 (10.26) 可得一个仅与速度 $S_L$ 与 $S_R$ 有关的 $S_*$ 公式：

$$
S_* = \frac{p_R - p_L + \rho_L u_L(S_L - u_L) - \rho_R u_R(S_R - u_R)}{\rho_L(S_L - u_L) - \rho_R(S_R - u_R)}.
\tag{10.27}
$$

于是，与较简单的 HLL 求解器一样，只需给出 $S_L$ 与 $S_R$ 的估计即可。

对 (10.23)、(10.24) 作代数运算并利用 (10.26) 的相应值 $p_L^*$、$p_R^*$，得中间通量 $\mathbf{F}_L^*$ 与 $\mathbf{F}_R^*$ 为

$$
\mathbf{F}_K^* = \mathbf{F}_K + S_K(\mathbf{U}_K^* - \mathbf{U}_K),
\tag{10.28}
$$

（$K = L$ 与 $K = R$），其中中间状态为

$$
\mathbf{U}_K^* = \rho_K\frac{S_K - u_K}{S_K - S_*}
\begin{bmatrix}
1\\[0.5ex]
S_*\\[0.5ex]
v_K\\[0.5ex]
w_K\\[1ex]
\dfrac{E_K}{\rho_K} + (S_* - u_K)\left[\dfrac{p_K}{\rho_K(S_K - u_K)} + S_*\right]
\end{bmatrix}.
\tag{10.29}
$$

HLLC 通量的最终选择按 (10.18) 作出。

下面介绍 HLLC 黎曼求解器的一个变体。由 (10.23) 与 (10.24) 可写出状态向量 $\mathbf{U}_L^*$ 与 $\mathbf{U}_R^*$ 的解：

$$
\mathbf{D}^* = [0, 1, 0, 0, S_*]^T,\qquad
\mathbf{U}_K^* = \frac{S_K\mathbf{U}_K - \mathbf{F}_K + p_K^*\mathbf{D}^*}{S_K - S_*},
\tag{10.30}
$$

其中 $p_L^*$ 与 $p_R^*$ 由 (10.26) 给出。把 (10.26) 的 $p_K^*$ 代入 (10.30)，再利用 (10.19) 与 (10.21)，便得中间通量的直接表达式：

$$
\mathbf{F}_K^* = \frac{S_*(S_K\mathbf{U}_K - \mathbf{F}_K) + S_K\big(p_K + \rho_K(S_K - u_K)(S_* - u_K)\big)\mathbf{D}^*}{S_K - S_*},
\tag{10.31}
$$

HLLC 通量的最终选择仍按 (10.18) 作出。

这里我们指出：HLLC 表述 (10.28)–(10.29) 强加了 $p_L^* = p_R^*$ 这一精确解所应满足的条件；而在另一种表述 (10.31) 中，我们放松了该条件，从而更符合压力近似 (10.26)。

若在星区内设单一的平均压力值，取 (10.26) 两压力的算术平均，即

$$
P_{LR} = \tfrac{1}{2}\left[p_L + p_R + \rho_L(S_L - u_L)(S_* - u_L) + \rho_R(S_R - u_R)(S_* - u_R)\right],
\tag{10.32}
$$

则得又一种 HLLC 通量。此时中间状态向量为

$$
\mathbf{U}_K^* = \frac{S_K\mathbf{U}_K - \mathbf{F}_K + P_{LR}\mathbf{D}^*}{S_K - S_*},
\tag{10.33}
$$

把它们代入 (10.19) 与 (10.21)，得通量 $\mathbf{F}_L^*$ 与 $\mathbf{F}_R^*$ 为

$$
\mathbf{F}_K^* = \frac{S_*(S_K\mathbf{U}_K - \mathbf{F}_K) + S_KP_{LR}\mathbf{D}^*}{S_K - S_*};
\tag{10.34}
$$

HLLC 通量的最终选择仍按 (10.18) 作出。

> **注** 迄今全部推导（只要波速 $S_L$ 与 $S_R$ 的估计可得）对**任意状态方程**均成立；状态方程仅进入 $S_L$ 与 $S_R$ 估计的规定之中。
{: .prompt-tip }

### 10.4.3 多维流动与多组分流动（Multidimensional and Multicomponent Flow）

这里考虑把 HLLC 求解器推广到两类应用：多维流动与多组分流动。

前文对 HLLC 格式的介绍是针对 x 方向分裂的三维 Euler 方程作出的，其中 $u$ 为法向速度。在一般的多维情形，界面未必与任何笛卡尔坐标方向对齐，但仍有一个法向速度分量与两个切向速度分量，前文所得全部结果均适用。

对多组分流动，组分可以与连续性方程联立写出守恒方程
$$
(\rho q_l)_t + (\rho u q_l)_x = 0,\qquad l = 1, \ldots, m .
$$

对扩充系统，其特征值不变，只是 $\lambda_2 = u$ 在三个空间维数下现在的重数为 $m + 3$。HLLC 通量无须任何特殊处理，被动标量可以和切向速度 $v$、$w$ 进行相同处理。

## 10.5 波速估计（Wave–Speed Estimates）

### 10.5.1 直接波速估计（Direct Wave Speed Estimates）

**1. 简单估计**

$$
S_L = u_L - a_L,\qquad S_R = u_R + a_R
\tag{10.35}
$$

$$
S_L = \min\{u_L - a_L,\ u_R - a_R\},\qquad S_R = \max\{u_L + a_L,\ u_R + a_R\}.
\tag{10.36}
$$

这些直接估计不推荐用于实际计算。

**2. Roe 平均特征值估计**

$$
S_L = \tilde{u} - \tilde{a},\qquad S_R = \tilde{u} + \tilde{a},
\tag{10.37}
$$

其中 $\tilde{u}$ 与 $\tilde{a}$ 分别为 Roe 平均的质点速度与声速，由下式给出：

$$
\tilde{u} = \frac{\sqrt{\rho_L}u_L + \sqrt{\rho_R}u_R}{\sqrt{\rho_L} + \sqrt{\rho_R}},\qquad
\tilde{a} = \left[(\gamma-1)\left(\tilde{H} - \tfrac{1}{2}\tilde{u}^2\right)\right]^{1/2},
\tag{10.38}
$$

其中焓 $H = (E+p)/\rho$ 近似为

$$
\tilde{H} = \frac{\sqrt{\rho_L}H_L + \sqrt{\rho_R}H_R}{\sqrt{\rho_L} + \sqrt{\rho_R}}.
\tag{10.39}
$$

**3. Einfeldt 波速**

在 Roe 特征值的启发下，Einfeldt 为其 HLLE 求解器提出了估计

$$
S_L = \bar{u} - \bar{d},\qquad S_R = \bar{u} + \bar{d},
\tag{10.40}
$$

其中

$$
\bar{d}^2 = \frac{\sqrt{\rho_L}a_L^2 + \sqrt{\rho_R}a_R^2}{\sqrt{\rho_L} + \sqrt{\rho_R}} + \eta^2(u_R - u_L)^2
\tag{10.41}
$$

且

$$
\eta^2 = \frac{1}{2}\,\frac{\sqrt{\rho_L}\sqrt{\rho_R}}{\left(\sqrt{\rho_L} + \sqrt{\rho_R}\right)^2}.
\tag{10.42}
$$

**4. Davis 观察**

对给定的黎曼问题能够确定一个正速度 $S^+$。在 HLL 通量中取 $S_L = -S^+$、$S_R = S^+$，便得到 Rusanov 通量：

$$
\mathbf{F}_{i+\frac{1}{2}} = \tfrac{1}{2}(\mathbf{F}_L + \mathbf{F}_R) - \tfrac{1}{2}S^+(\mathbf{U}_R - \mathbf{U}_L).
$$

当

$$
S^+ = \max\{|u_L| + a_L,\ |u_R| + a_R\};
$$

鲁棒性强。若取

$$
S^+=S^n_{max} = \frac{C_{cfl}\Delta x}{\Delta t}
$$

当 $C_{cfl} = 1$ 时 $S^+ = \Delta x/\Delta t$，相应得到 Lax–Friedrichs 数值通量

$$
\mathbf{F}_{i+\frac{1}{2}} = \tfrac{1}{2}(\mathbf{F}_L + \mathbf{F}_R) - \tfrac{1}{2}\frac{\Delta x}{\Delta t}(\mathbf{U}_R - \mathbf{U}_L).
$$

### 10.5.2 基于压力的波速估计（Pressure–Based Wave Speed Estimates）

选取波速

$$
S_L = u_L - a_L q_L,\qquad S_R = u_R + a_R q_R,
\tag{10.43}
$$

其中

$$
q_K = \begin{cases}
1 & \text{if } p^* \le p_K \text{（稀疏波）},\\[1ex]
\left[1 + \dfrac{\gamma+1}{2\gamma}\left(\dfrac{p^*}{p_K} - 1\right)\right]^{1/2} & \text{if } p^* > p_K \text{（激波）} .
\end{cases}
\tag{10.44}
$$

这种波速选取能够区分激波与稀疏波。建议采用第 9 章的压力近似来求 $p^*$。第 9 章第 9.3 节给出的 PVRS 近似黎曼求解器给出

$$
p_{pvrs} = \tfrac{1}{2}(p_L + p_R) - \tfrac{1}{2}(u_R - u_L)\bar{\rho}\bar{a},
\tag{10.45}
$$

其中

$$
\bar{\rho} = \tfrac{1}{2}(\rho_L + \rho_R),\qquad \bar{a} = \tfrac{1}{2}(a_L + a_R).
\tag{10.46}
$$

$p^*$ 的另一选择由第 9 章第 9.4.1 节的双稀疏波黎曼求解器 TRRS 提供，即

$$
p_{tr} = \left\{\frac{a_L + a_R - \dfrac{\gamma-1}{2}(u_R - u_L)}{\dfrac{a_L}{p_L^{z}} + \dfrac{a_R}{p_R^{z}}}\right\}^{1/z},
\tag{10.63}
$$

其中

$$
z = \frac{\gamma-1}{2\gamma},\qquad P_{LR} = \frac{p_L}{p_R}.
\tag{10.64}
$$

> 原书在第 10 章前文中将符号 $P_{LR}$ 用作星区平均压力，而此处沿用第 9 章的记号将 $P_{LR}$ 定义为压力比 $p_L/p_R$，二者为符号重用，请读者注意区分。
{: .prompt-tip }

第 9 章第 9.4.2 节的双激波黎曼求解器 TSRS 给出

$$
p_{ts} = \frac{g_L(p_0)p_L + g_R(p_0)p_R - \Delta u}{g_L(p_0) + g_R(p_0)},
\tag{10.65}
$$

其中

$$
g_K(p) = \left(\frac{A_K}{p + B_K}\right)^{1/2},\qquad p_0 = \max(0, p_{pvrs}),
\tag{10.66}
$$

## 10.6 HLLC 通量小结（Summary of HLLC Fluxes）

这里基于一种特定的波速选取总结 HLLC 格式。计算 HLLC 通量按如下步骤进行。

**步骤 I：压力估计。** 按下式计算星区压力 $p^*$ 的估计：

$$
p^* = \max(0, p_{pvrs}),\qquad
p_{pvrs} = \tfrac{1}{2}(p_L + p_R) - \tfrac{1}{2}(u_R - u_L)\bar{\rho}\bar{a},
$$

其中 $\bar{\rho} = \tfrac{1}{2}(\rho_L + \rho_R)$，$\bar{a} = \tfrac{1}{2}(a_L + a_R)$。

**步骤 II：波速估计。** 按下式计算 $S_L$ 与 $S_R$ 的估计：

$$
S_L = u_L - a_L q_L,\qquad S_R = u_R + a_R q_R,
$$

其中

$$
q_K = \begin{cases}
1 & \text{if } p^* \le p_K,\\[1ex]
\left[1 + \dfrac{\gamma+1}{2\gamma}\left(\dfrac{p^*}{p_K} - 1\right)\right]^{1/2} & \text{if } p^* > p_K ;
\end{cases}
$$

再按下式由 $S_L$ 与 $S_R$ 计算中间速度 $S_*$：

$$
S_* = \frac{p_R - p_L + \rho_L u_L(S_L - u_L) - \rho_R u_R(S_R - u_R)}{\rho_L(S_L - u_L) - \rho_R(S_R - u_R)}.
$$

**步骤 III：HLLC 通量。** 按下式计算 HLLC 通量：

$$
\mathbf{F}_{i+\frac{1}{2}}^{hllc} = \begin{cases}
\mathbf{F}_L & \text{if } 0 \le S_L,\\
\mathbf{F}_L^* & \text{if } S_L \le 0 \le S_*,\\
\mathbf{F}_R^* & \text{if } S_* \le 0 \le S_R,\\
\mathbf{F}_R & \text{if } 0 \ge S_R ,
\end{cases}
$$

其中

$$
\mathbf{F}_K^* = \mathbf{F}_K + S_K(\mathbf{U}_K^* - \mathbf{U}_K)
$$

且

$$
\mathbf{U}_K^* = \rho_K\frac{S_K - u_K}{S_K - S_*}
\begin{bmatrix}
1\\[0.5ex]
S_*\\[0.5ex]
v_K\\[0.5ex]
w_K\\[1ex]
\dfrac{E_K}{\rho_K} + (S_* - u_K)\left[\dfrac{p_K}{\rho_K(S_K - u_K)} + S_*\right]
\end{bmatrix}.
$$

如上所见，第三步中的 HLLC 通量有两种变体。

**步骤 III：HLLC 通量，变体 1。** 按下式计算数值通量：

$$
\mathbf{F}_K^* = \frac{S_*(S_K\mathbf{U}_K - \mathbf{F}_K) + S_K\big(p_K + \rho_K(S_K - u_K)(S_* - u_K)\big)\mathbf{D}^*}{S_K - S_*},\\
\mathbf{D}^* = [0, 1, 0, 0, S_*]^T,
$$

**步骤 III：HLLC 通量，变体 2。** 按下式计算数值通量：

$$
\mathbf{F}_K^* = \frac{S_*(S_K\mathbf{U}_K - \mathbf{F}_K) + S_KP_{LR}\mathbf{D}^*}{S_K - S_*},
$$

其中 $\mathbf{D}^*$ 同上，而

$$
P_{LR} = \tfrac{1}{2}\left[p_L + p_R + \rho_L(S_L - u_L)(S_* - u_L) + \rho_R(S_R - u_R)(S_* - u_R)\right].
$$

## 10.7 接触波与被动标量（Contact Waves and Passive Scalars）

本节从略，其要点为：HLL 格式在分辨物质界面、剪切波与涡等特殊但重要的流动特征时，会引入过量的数值耗散。

## 10.8 数值结果（Numerical Results）

> 本节从略。原书在一维与多维标准算例上对 HLL 与 HLLC 的数值表现进行了对比。
{: .prompt-info }
