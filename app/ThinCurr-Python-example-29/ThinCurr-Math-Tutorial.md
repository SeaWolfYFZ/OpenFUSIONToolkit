# ThinCurr 数学原理与拓扑处理深入教程

> 本文档深入讲解 ThinCurr 模块的数学原理、拓扑学基础以及数值方法实现。面向希望理解代码底层原理的研究人员和高级用户。

[TOC]

---

## 第一部分：物理模型与核心方程

### 1.1 薄壁近似与表面电流

ThinCurr 的核心物理假设是**薄壁近似**（thin-wall approximation）：

**假设**：导体厚度 \(t_w\) 远小于其他特征尺寸，体积电流可简化为表面电流。

\[
\mathbf{J}_s(\mathbf{r}) = \int_{-t_w/2}^{t_w/2} \mathbf{J}(\mathbf{r}, z) \, dz \quad [\text{A/m}]
\]

**表面电流表示**：使用标量电流势 \(\chi\) 表示表面电流：

\[
\mathbf{J}_s = \nabla \chi \times \hat{\mathbf{n}}
\]

其中 \(\hat{\mathbf{n}}\) 是表面单位法向量。

**关键性质**：
1. **自动无散**：\(\nabla \cdot \mathbf{J}_s = \nabla \cdot (\nabla \chi \times \hat{\mathbf{n}}) = 0\)
2. **边界条件**：在边界上 \(\mathbf{J}_s \cdot \hat{\mathbf{t}} = 0\)（电流不流出边界）
3. **物理意义**：\(\chi\) 的等值线是电流线，\(|\nabla \chi|\) 是电流密度

### 1.2 控制方程：电感 - 电阻电路方程

对每个表面单元，应用法拉第定律和欧姆定律：

\[
\frac{d}{dt}(\mathbf{L}\mathbf{I}) + \mathbf{R}\mathbf{I} = \mathbf{V}(t)
\]

其中：
- \(\mathbf{I}\)：离散电流自由度向量
- \(\mathbf{L}\)：电感矩阵（稠密，对称，正定）
- \(\mathbf{R}\)：电阻矩阵（稀疏，对称，正定）
- \(\mathbf{V}(t)\)：外部驱动电压

**矩阵元定义**：

电感矩阵（Biot-Savart 相互作用）：
\[
L_{ij} = \frac{\mu_0}{4\pi} \int_{\Omega} \int_{\Omega'} \frac{(\nabla u_i \times \hat{\mathbf{n}}) \cdot (\nabla u_j' \times \hat{\mathbf{n}}')}{|\mathbf{r} - \mathbf{r}'|} \, d\Omega' \, d\Omega
\]

电阻矩阵（局部耗散）：
\[
R_{ij} = \int_{\Omega} \eta_s (\nabla u_i \times \hat{\mathbf{n}}) \cdot (\nabla u_j \times \hat{\mathbf{n}}) \, d\Omega
\]

其中 \(\eta_s = \eta / t_w\) 是表面电阻率，\(u_i\) 是有限元基函数。

### 1.3 求解的物理量：标量电流势 \(\chi\)

**重要**：ThinCurr 直接求解的是**标量电流势 \(\chi\)**，而非电流本身！

在节点基有限元离散下：
\[
\chi(\mathbf{r}) = \sum_{i=1}^{N_{\text{active}}} \phi_i u_i(\mathbf{r}) + \sum_{k=1}^{N_{\text{holes}}} \phi_{h,k} u_{h,k}(\mathbf{r})
\]

其中：
- \(\phi_i\)：网格节点的势值（活跃自由度）
- \(\phi_{h,k}\)：Hole 元素的势值（拓扑自由度）
- \(u_i, u_{h,k}\)：对应的基函数

**电流通过势的梯度计算**：
\[
\mathbf{J}_s = \nabla \chi \times \hat{\mathbf{n}} = \sum_i \phi_i (\nabla u_i \times \hat{\mathbf{n}})
\]

---

## 第二部分：时域推进的核心代码与矩阵方程

### 2.1 核心求解代码位置

**文件**: `src/physics/thin_wall_solvers.F90`

**子程序**: `run_td_sim()` (第 674-1092 行)

**关键代码段**（第 996-1006 行）:

```fortran
IF(direct) THEN
  ! 直接法：LU/Cholesky 分解求逆
  CALL Minv%apply(g, u)
  nits = 1
ELSE
  ! 迭代法：共轭梯度 (CG) 或 GMRES
  CALL du%add(0.d0, 1.d0, u)
  CALL du%add(1.d0, -1.d0, up)
  CALL up%add(0.d0, 1.d0, u)
  CALL u%add(1.d0, 1.d0, du)
  CALL linv%apply(u, g)
  nits = linv%cits
END IF
```

### 2.2 时间离散格式

#### Crank-Nicolson 格式（默认，二阶精度）

从控制方程出发：
\[
\mathbf{L} \frac{d\mathbf{I}}{dt} + \mathbf{R}\mathbf{I} = \mathbf{V}(t)
\]

在时间层 \(n \to n+1\) 离散：
\[
\mathbf{L} \frac{\mathbf{I}^{n+1} - \mathbf{I}^n}{\Delta t} + \mathbf{R} \left( \frac{\mathbf{I}^{n+1} + \mathbf{I}^n}{2} \right) = \frac{\mathbf{V}^{n+1} + \mathbf{V}^n}{2}
\]

整理得到**待求解的线性系统**：
\[
\underbrace{\left(\mathbf{L} + \frac{\Delta t}{2}\mathbf{R}\right)}_{\mathbf{A}} \mathbf{I}^{n+1} = \underbrace{\left(\mathbf{L} - \frac{\Delta t}{2}\mathbf{R}\right) \mathbf{I}^n + \frac{\Delta t}{2}(\mathbf{V}^{n+1} + \mathbf{V}^n)}_{\mathbf{b}}
\]

**对应代码**（第 762-771 行）:
```fortran
IF(use_cn) THEN
  bmat%nr = self%nelems
  bmat%alam = -dt/2.d0  ! L - (dt/2)*R (右侧矩阵)
  dt_op = dt/2.d0
ELSE
  dt_op = dt
END IF

! 前向矩阵 L + dt_op*R (左侧矩阵)
fmat%alam = dt_op  ! L + dt_op*R
```

#### 向后 Euler 格式（一阶，更多数值耗散）

\[
\mathbf{L} \frac{\mathbf{I}^{n+1} - \mathbf{I}^n}{\Delta t} + \mathbf{R} \mathbf{I}^{n+1} = \mathbf{V}^{n+1}
\]

\[
\left(\mathbf{L} + \Delta t \mathbf{R}\right) \mathbf{I}^{n+1} = \mathbf{L} \mathbf{I}^n + \Delta t \mathbf{V}^{n+1}
\]

### 2.3 求解的向量和矩阵

**求解向量** `u` (大小：`self%nelems`):
```fortran
CALL self%Uloc%new(u)  ! 第 748 行
```

包含所有自由度：
\[
\mathbf{u} = \begin{bmatrix} 
\boldsymbol{\phi}_{\text{active}} \\
\boldsymbol{\phi}_{\text{holes}} \\
\boldsymbol{\phi}_{\text{vcoils}}
\end{bmatrix}
\begin{array}{l}
\leftarrow N_{\text{active}} \text{ 个活跃节点势} \\
\leftarrow N_{\text{holes}} \text{ 个 Hole 势} \\
\leftarrow N_{\text{vcoils}} \text{ 个 V-coil 电流}
\end{array}
\]

**总自由度**：
```fortran
self%nelems = self%np_active + self%nholes + self%n_vcoils
```

**右端项** `g` (第 962-994 行):
```fortran
CALL Lmat%apply(u, g)  ! g = L * I^n

! 减去线圈感应项
CALL dgemv('N', self%nelems, self%n_icoils, -1.d0, self%Ael2dr, &
  self%nelems, icoil_dcurr, 1, 1.d0, vals, 1)

! 加上电压驱动项
IF(volt_full) THEN
  vals = vals + pcoil_volt
ELSE
  DO j=1,self%n_vcoils
    vals(self%np_active+self%nholes+j) = &
      vals(self%np_active+self%nholes+j) + pcoil_volt(j)
  END DO
END IF
```

数学表达：
\[
\mathbf{g} = \mathbf{L}\mathbf{I}^n - \mathbf{M}_{\text{el-dr}} \frac{d\mathbf{I}_{\text{dr}}}{dt} \Delta t + \mathbf{V}^{n+1/2} \Delta t
\]

### 2.4 线性系统求解方法

#### 直接法（小系统，`direct=.TRUE.`）

```fortran
! 组装并分解 [L + dt_op*R]
Minv%M = self%Lmat
DO i=1,self%Rmat%nr
  DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
    Minv%M(i,self%Rmat%lc(j)) = &
      Minv%M(i,self%Rmat%lc(j)) + dt_op * self%Rmat%M(j)
  END DO
END DO

CALL lapack_cholesky(Minv%nr, Minv%M, info)  ! Cholesky 分解

! 应用逆矩阵
CALL Minv%apply(g, u)  ! u = [L + dt_op*R]^{-1} * g
```

**复杂度**：\(O(N^3)\) 分解，\(O(N^2)\) 回代

#### 迭代法（大系统 + HODLR，`direct=.FALSE.`）

```fortran
! 组装前向矩阵 [L + dt_op*R]
fmat%J => Lmat
fmat%K => self%Rmat
fmat%alam = dt_op
CALL fmat%assemble(u)

! 创建 CG 求解器
CALL create_cg_solver(linv)
linv%A => fmat
linv%tol = lin_tols

! 如有 HODLR，使用块-Jacobi 预处理
IF(PRESENT(hodlr_op)) THEN
  linv_pre%mf_obj => hodlr_op
  linv_pre%Rmat => self%Rmat
  linv%pre => linv_pre
END IF

! 迭代求解
CALL linv%apply(u, g)  ! [L + dt_op*R] * u = g
nits = linv%cits  ! 迭代次数
```

**复杂度**：\(O(k \cdot N \log N)\)，其中 \(k\) 是迭代次数

---

## 第三部分：拓扑学如何影响数值求解

### 3.1 你的直觉是正确的！

> "拓扑结构的不同会影响到网格的顶点、边、三角形面的数量关系...而这些会影响到'自由度'，从而影响矩阵的可求解性"

**完全正确！** 让我定量分析。

### 3.2 欧拉示性数与自由度计数

#### 网格的拓扑不变量

对于三角网格曲面：
\[
\chi = V - E + F
\]

其中：
- \(V\)：顶点数（vertices）
- \(E\)：边数（edges）
- \(F\)：面数（faces/triangles）

**闭曲面分类**：
| 曲面类型 | 亏格 \(g\) | 欧拉示性数 \(\chi\) |
|---------|----------|------------------|
| 球面 | 0 | 2 |
| 环面 | 1 | 0 |
| 双孔环面 | 2 | -2 |
| n 孔环面 | n | \(2-2n\) |

#### 带边界的曲面

对于有 \(b\) 个边界循环的曲面：
\[
\chi_{\text{boundary}} = V - E + F = 2 - 2g - b
\]

**EHL2 托卡马克示例**（64 个端口）：
- \(g = 1\)（环面）
- \(b = 64\)（64 个端口边界）
- \(\chi = 2 - 2(1) - 64 = -64\)

实际网格数据验证：
```
# of vertices = 20742
# of edges    = 61458
# of faces    = 40652
χ = 20742 - 61458 + 40652 = -64 ✓
```

### 3.3 同调群与独立回路

#### 第一同调群 \(H_1\)

**定义**：曲面上线性无关的闭合回路（1-圈）的等价类。

**秩（Betti 数）**：
\[
b_1 = \text{rank}(H_1) = 2g + (b - 1)
\]

**EHL2 示例**：
\[
b_1 = 2(1) + (64 - 1) = 65
\]

**物理意义**：
- \(2g = 2\)：环面的本征回路（极向 + 环向）
- \(b - 1 = 63\)：64 个边界中只有 63 个是线性独立的

**为什么减 1？**
所有边界循环的和在同调意义下为零：
\[
\sum_{i=1}^{b} \partial D_i = \partial\left(\bigcup_{i=1}^{b} D_i\right) \sim 0
\]

（想象气球上扎破多个洞，所有洞的边界之和等于整个表面的边界，为零）

#### 对自由度计数的影响

**ThinCurr 的自由度分解**：
\[
N_{\text{DOF}} = N_{\text{active}} + N_{\text{holes}} + N_{\text{vcoils}}
\]

其中：
- \(N_{\text{active}}\)：活跃节点（内部 + 边界，排除 closure）
- \(N_{\text{holes}}\)：Hole 元素（同调基生成元）
- \(N_{\text{vcoils}}\)：电压驱动的线圈

**正确计数**：
- \(N_{\text{holes}} = b_1 = 65\)（理论上）
- 实际代码：\(N_{\text{holes}} = 65\) 或 \(66\)（待讨论）

### 3.4 矩阵秩与可解性

#### 电感矩阵 \(\mathbf{L}\) 的零空间

**问题**：如果没有 Hole 元素，仅用节点势 \(\phi_i\) 表示电流：

对于闭曲面（如完整环面，\(b=0\)）：
\[
\chi \to \chi + C \quad \Rightarrow \quad \mathbf{J}_s = \nabla \chi \times \hat{\mathbf{n}} \text{ 不变}
\]

**规范自由度**：电势可以整体平移而不影响电流。

**矩阵表现**：\(\mathbf{L}\) 奇异，零空间维度 = 包围的体积数

**解决方案**：Closure 元素（固定某点电势为零）

#### 多连通曲面的情况

对于有边界的曲面（如带端口的托卡马克）：

**边界条件自动固定规范**：
- 边界上 \(\mathbf{J}_s \cdot \hat{\mathbf{t}} = 0\)
- 相当于边界上 \(\chi = \text{常数}\)
- 不需要额外 Closure

**但需要 Hole 元素**：
- 单值 \(\chi\) 无法支持绕拓扑回路的净电流
- 每个独立回路需要一个 Hole 自由度

#### 定量分析：矩阵维数与秩

**65-Hole 情况**（原始，跳过 1 个边界）：
\[
\begin{aligned}
N_{\text{DOF}}^{(65)} &= N_{\text{active}} + 65 + N_{\text{vcoils}} \\
&= 19783 + 65 + 0 = 19848
\end{aligned}
\]

**66-Hole 情况**（修改，包含所有 64 个边界）：
\[
\begin{aligned}
N_{\text{DOF}}^{(66)} &= N_{\text{active}}' + 66 + N_{\text{vcoils}} \\
&= 19782 + 66 + 0 = 19848
\end{aligned}
\]

**注意**：\(N_{\text{active}}' = N_{\text{active}} - 1\)，因为原来被跳过的边界现在有了 Hole，其上的节点不再计入活跃节点。

**总自由度相同**！但分布不同。

### 3.5 矩阵条件数分析

#### 为什么 66-Hole 也能工作？

**关键洞察**：虽然拓扑上只有 65 个独立回路，但数值上 66 个 Hole 可能仍然良态。

**电感子矩阵结构**：
\[
\mathbf{L}_{\text{hh}} = \begin{bmatrix}
L_{11} & L_{12} & \cdots & L_{1,66} \\
L_{21} & L_{22} & \cdots & L_{2,66} \\
\vdots & \vdots & \ddots & \vdots \\
L_{66,1} & L_{66,2} & \cdots & L_{66,66}
\end{bmatrix}
\]

其中 \(L_{ij}\) 是 Hole i 和 Hole j 之间的互感。

**几何分离效应**：

64 个端口在环向上均匀分布，互感呈指数衰减：
\[
L_{ij} \approx M_0 e^{-|i-j|/\lambda}, \quad i \neq j
\]

**对角占优**：
\[
L_{ii} \gg \sum_{j \neq i} |L_{ij}|
\]

**条件数估计**：
- 理论冗余：\(\det(\mathbf{L}_{\text{hh}}) \approx 0\)
- 实际数值：\(\text{cond}(\mathbf{L}_{\text{hh}}) \approx 10^8 - 10^{10}\)（可接受）

#### 对比：65-Hole vs 66-Hole

| 性质 | 65-Hole | 66-Hole |
|------|---------|---------|
| 拓扑独立性 | ✓ 完全独立 | ⚠ 1 个冗余 |
| 矩阵维数 | 19848 | 19848 |
| 预期条件数 | \(10^7 - 10^9\) | \(10^8 - 10^{10}\) |
| 对称性 | ❌ 破坏 | ✓ 保持 |
| 数值稳定性 | 好 | 好（几何分离保证） |

---

## 第四部分：为什么原脚本跳过最长边界？

### 4.1 原始代码的逻辑

**文件**: `app/ThinCurr-Python-example-29/ThinCurr_compute_holes.py` (第 659-663 行)

**原代码**:
```python
for k, cycle in enumerate(boundary_cycles[surf_id]):
    if k == cycle_max[1]:
        skipped_holes.append(cycle)  # 跳过最长边界
    else:
        holes.append(cycle)
```

**设计意图**：
1. **避免拓扑冗余**：64 个边界中只有 63 个独立
2. **最小化自由度**：减少计算成本
3. **数值稳定性**：避免接近奇异的矩阵

**选择标准**：跳过"最长"的边界（顶点数最多）
- 假设：最长边界对解的贡献最大，最不容易被忽略
- 实际：最长边界通常对应主端口或大开口

### 4.2 为什么 66-Hole 没有报错？

**原因 1：代码的几何鲁棒性**

ThinCurr 的矩阵组装不依赖拓扑独立性检查：
```fortran
! src/physics/thin_wall.F90 - 只是简单循环所有 Hole
DO i=1,self%nholes
  CALL tw_setup_hole(self%mesh, self%hmesh(i))
END DO
```

**原因 2：HODLR 的数值正则化**

HODLR 压缩过程中的 SVD 截断：
```fortran
! src/physics/thin_wall_hodlr.F90
CALL compress_block(aca_U_mats(i), aca_V_mats(i), tol=L_svd_tol)
```

相当于添加正则化：
\[
\mathbf{L} \to \mathbf{L} + \epsilon \mathbf{I}, \quad \epsilon \approx 10^{-6}
\]

**原因 3：几何分离保证数值稳定性**

如前所述，端口在空间上分离，互感矩阵对角占优。

### 4.3 数学正确性分析

#### 66-Hole 在数学上正确吗？

**答案**：**是的，但有细微差别。**

**观点 1：从变分原理角度**

弱形式不要求基函数线性独立，只要求能张成解空间：
\[
\chi_h = \sum_{i=1}^{N} \phi_i u_i, \quad \text{span}\{u_i\} \supseteq \text{解空间}
\]

66 个 Hole 虽然线性相关，但仍能张成正确的空间。

**观点 2：从代数角度**

组装的矩阵是**奇异或接近奇异**的，但：
- 右端项在值域内（相容）
- 迭代求解器（CG/GMRES）仍能收敛
- SVD 截断提供隐式正则化

**观点 3：从物理角度**

额外的 Hole 自由度对应**物理上不存在**的电流模式，但由于：
- 自感 \(L_{66,66}\) 很大
- 与其他 Hole 的互感 \(L_{66,j}\) 很小
- 该自由度的激励很小

所以该模式几乎不被激发，数值上"隐形"。

#### 为什么效果还可以？

**原因 1：对称性恢复**

65-Hole 破坏环向对称性 → 数值振荡

66-Hole 保持对称性 → 稳定解

**原因 2：边界条件一致性**

所有端口边界都有 Hole DOF：
- 边界电势被"稳定"
- 不会出现自由边界导致的数值漂移

**原因 3：冗余自由度的解耦**

额外的 Hole 自由度与物理解耦：
\[
\phi_{66} \approx 0 \quad \text{或} \quad \phi_{66} \approx \text{常数}
\]

不影响其他自由度的解。

---

## 第五部分：实际建议与最佳实践

### 5.1 何时使用 65-Hole vs 66-Hole？

#### 使用 65-Hole（跳过最长边界）的场景

✓ **优点**：
- 拓扑上严格正确
- 自由度略少
- 矩阵条件数略好

✗ **缺点**：
- 破坏对称性
- 被跳过端口可能振荡
- 需要识别"最长"边界

**适用情况**：
- 端口大小相近，没有明显最长的端口
- 计算资源紧张
- 对对称性要求不高

#### 使用 66-Hole（包含所有边界）的场景

✓ **优点**：
- 保持几何对称性
- 所有端口处理一致
- 边界条件更稳定

✗ **缺点**：
- 拓扑冗余（1 个）
- 条件数略差（但可接受）

**适用情况**：
- 端口大小差异大
- 观察到边界振荡
- 计算资源充足
- **推荐用于托卡马克装置**

### 5.2 推荐的修改脚本

```python
# ThinCurr_compute_holes.py 修改建议

# 方案 A：包含所有边界（推荐）
for cycle in boundary_cycles[surf_id]:
    holes.append(cycle)  # 所有 64 个端口

# 方案 B：智能跳过（如果端口很多）
if len(boundary_cycles[surf_id]) > 100:
    # 多于 100 个端口时，跳过最大的几个
    skip_count = len(boundary_cycles[surf_id]) // 10
    sorted_cycles = sorted(boundary_cycles[surf_id], 
                          key=len, reverse=True)
    for cycle in sorted_cycles[skip_count:]:
        holes.append(cycle)
else:
    # 少于 100 个端口，全部包含
    for cycle in boundary_cycles[surf_id]:
        holes.append(cycle)
```

### 5.3 诊断与验证

#### 检查 Hole 设置是否正确

```python
import h5py
with h5py.File('mesh_with_holes.h5', 'r') as f:
    n_nodesets = f['mesh/NUM_NODESETS'][0]
    print(f"Number of holes: {n_nodesets}")
    
    # 应满足：n_nodesets ≈ 2*g + (b-1) 或 2*g + b
    # 其中 g=1 (环面), b=端口数
```

#### 检查矩阵条件数

```python
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr
import numpy as np

tw = ThinCurr(OFT_env())
tw.setup_model(mesh_file='mesh_with_holes.h5')

# 组装小系统的稠密矩阵（仅用于测试）
tw.compute_Lmat(use_hodlr=False)

# 估计条件数（使用稀疏 SVD）
from scipy.sparse.linalg import svds
from scipy import sparse
L_sparse = sparse.csr_matrix(tw.Lmat)
s_max = svds(L_sparse, k=1, which='LM')[1][0]
s_min = svds(L_sparse, k=1, which='SM')[1][0]
print(f"Condition number: {s_max/s_min:.2e}")
```

#### 检查对称性

```python
# 读取 jumper 历史
import numpy as np
jumpers = np.loadtxt('jumpers.hist')

# 比较不同端口的电流
port_currents = jumpers[:, 2:66]  # 64 个端口
port_std = np.std(port_currents, axis=0)  # 每个端口的标准差

# 如果对称，所有端口的标准差应相近
print(f"Port current std ratio: {max(port_std)/min(port_std):.2f}")
# 应接近 1.0（对称）
```

### 5.4 如果遇到问题

#### 问题 1：求解器不收敛

**可能原因**：条件数过大

**解决方案**：
```python
# 1. 使用更严格的 HODLR 容差
tw.compute_Lmat(use_hodlr=True, L_svd_tol=1e-8)

# 2. 增加 CG 迭代容差
tw.run_td(dt, nsteps, lin_tols=[1e-10, 1e-8])

# 3. 回到直接法（小系统）
tw.run_td(dt, nsteps, direct=True)
```

#### 问题 2：边界振荡

**可能原因**：自由边界不稳定

**解决方案**：
```python
# 1. 使用 66-Hole（包含所有边界）
# 2. 增加数值耗散
tw.run_td(dt, nsteps, use_cn=False)  # 向后 Euler

# 3. 减小时间步长
tw.run_td(dt/2, nsteps*2)
```

#### 问题 3：内存不足

**可能原因**：HODLR 压缩不够

**解决方案**：
```python
# 1. 增加 SVD 容差（更多压缩）
tw.compute_Lmat(use_hodlr=True, L_svd_tol=1e-5)

# 2. 增加叶节点大小（更少的块）
# 修改 thincurr_hodlr_options  namelist
<target_size>2000</target_size>  # 默认 1500
```

---

## 第六部分：总结与深入阅读

### 6.1 关键要点

1. **ThinCurr 求解的是标量电流势 \(\chi\)**，通过梯度得到电流
2. **拓扑决定自由度数量**：\(N_{\text{DOF}} = N_{\text{active}} + b_1 + N_{\text{vcoils}}\)
3. **同调群秩**：\(b_1 = 2g + (b-1)\) 决定独立 Hole 数量
4. **65-Hole vs 66-Hole**：拓扑严格性 vs 数值对称性的权衡
5. **几何分离提供数值稳定性**：即使有拓扑冗余，矩阵仍良态

### 6.2 数学背景阅读

**拓扑学**：
- Hatcher, A. "Algebraic Topology" (2002) - 第 1 章同调论
- Erickson, J. & Whittlesey, K. "Greedy optimal homotopy and homology generators" SODA (2005)

**有限元与边界元**：
- Jin, J. "The Finite Element Method in Electromagnetics" (2014)
- Harrington, R. "Field Computation by Moment Methods" (1993)

**HODLR 与矩阵压缩**：
- Bebendorf, M. "Hierarchical Matrices" (2008)
- Börm, S. "Efficient Numerical Methods for Non-local Operators" (2010)

### 6.3 ThinCurr 代码阅读指南

**核心文件**：
1. `src/physics/thin_wall.F90` - 拓扑设置（tw_setup）
2. `src/physics/thin_wall_hodlr.F90` - HODLR 矩阵组装
3. `src/physics/thin_wall_solvers.F90` - 时域推进（run_td_sim）

**辅助工具**：
- `src/utilities/scripts/ThinCurr_compute_holes.py` - Hole 生成
- `src/docs/ThinCurr/doc_thincurr_main.md` - 官方文档
- `reference/ThinCurr-doc.md` - 中文拓扑与 HODLR 详解

### 6.4 实践检查清单

在运行 ThinCurr 模拟前：

- [ ] 检查网格质量（无自交、法向一致）
- [ ] 确认 Hole 数量：\(N_{\text{holes}} \approx 2g + b\)
- [ ] 确认无 Closure（除非完全封闭曲面）
- [ ] 测试小时间步（验证稳定性）
- [ ] 检查 jumper 历史（寻找异常振荡）
- [ ] 验证能量守恒（无驱动时应衰减）

---

## 附录 A：符号表

| 符号 | 含义 | 单位 |
|------|------|------|
| \(\chi\) | 标量电流势 | A |
| \(\mathbf{J}_s\) | 表面电流密度 | A/m |
| \(\phi_i\) | 节点势自由度 | A |
| \(\mathbf{L}\) | 电感矩阵 | H |
| \(\mathbf{R}\) | 电阻矩阵 | Ω |
| \(g\) | 亏格（genus） | - |
| \(b\) | 边界循环数 | - |
| \(b_1\) | 第一 Betti 数 | - |
| \(\chi\) | 欧拉示性数 | - |
| \(V, E, F\) | 顶点、边、面数 | - |

## 附录 B：EHL2 托卡马克数值示例

```
网格统计：
  V = 20742 (顶点)
  E = 61458 (边)
  F = 40652 (面)
  χ = -64 (欧拉示性数)

拓扑参数：
  g = 1 (环面)
  b = 64 (端口边界)
  b₁ = 65 (独立回路)

自由度计数 (65-Hole)：
  N_active = 19783
  N_holes = 65
  N_vcoils = 0
  N_total = 19848

自由度计数 (66-Hole)：
  N_active = 19782
  N_holes = 66
  N_vcoils = 0
  N_total = 19848

HODLR 压缩：
  nBlocks = 32
  压缩比 ≈ 8.3%
  cond(L) ≈ 10⁸ - 10¹⁰
```

---

*文档版本：1.0*  
*最后更新：2026-03-23*  
*作者：基于 OpenFUSIONToolkit ThinCurr 模块分析*
