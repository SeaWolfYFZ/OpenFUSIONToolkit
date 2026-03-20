# ThinCurr 技术文档：网格拓扑处理与HODLR矩阵压缩

## 目录

1. [概述](#概述)
2. [网格拓扑处理](#网格拓扑处理)
   - 2.1 [物理背景与数学基础](#物理背景与数学基础)
   - 2.2 ["洞"(Hole)元素](#洞hole元素)
   - 2.3 ["闭合"(Closure)元素](#闭合closure元素)
   - 2.4 [代码实现](#代码实现)
3. [HODLR矩阵压缩](#hodlr矩阵压缩)
   - 3.1 [算法原理](#算法原理)
   - 3.2 [空间划分](#空间划分)
   - 3.3 [自适应交叉近似(ACA+)](#自适应交叉近似aca)
   - 3.4 [代码实现](#代码实现-1)
4. [时域模拟中的处理](#时域模拟中的处理)
5. [参考资源](#参考资源)

---

## 概述

ThinCurr是一个开源的边界有限元方法(BFEM)代码，用于模拟三维导电结构（如聚变装置的真空容器）中的涡流动力学。本文档详细介绍ThinCurr处理网格拓扑和实现HODLR矩阵压缩的技术细节，帮助用户理解其底层原理。

---

## 网格拓扑处理

### 物理背景与数学基础

在薄壁模型中，表面电流密度 $\mathbf{J}_s$ 通过标量势函数 $\chi$ 表示：

$$\mathbf{J}_s = \nabla\chi \times \hat{\mathbf{n}}$$

其中 $\hat{\mathbf{n}}$ 是表面的一致取向单位法向量。这种表示自动满足电流的无散条件（$\nabla \cdot \mathbf{J}_s = 0$），因为电流仅依赖于 $\chi$ 的梯度。

然而，这种表示方法在处理**多连通区域**（multiply-connected domains）时会遇到拓扑学问题：

1. **边界问题**：开曲面（如圆柱面）具有边界，边界上的电势必须为常数
2. **闭合问题**：闭曲面（如球面、环面）没有边界，但电势必须固定参考点
3. **多值问题**：环面等几何体存在拓扑上不可收缩的回路，单值电势无法支持这些回路中的净电流

### "洞"(Hole)元素

#### 拓扑学概念

在拓扑学中，"洞"是指几何体中使其成为多连通区域的特征。例如：
- 圆柱面有一个"洞"（管状结构）
- 环面（torus）有两个"洞"：极向和环向

#### ThinCurr中的实现

为了解决单值电势的限制，ThinCurr引入了**"洞"元素**作为特殊的自由度。这些元素代表在拓扑不可收缩回路上施加的常数电势，使得：

1. **允许多值电势**：通过"洞"元素，电势 $\chi$ 可以在封闭回路上产生差值
2. **支持净电流**：可以模拟围绕"洞"边界流动的总电流
3. **捕获磁通**：捕获通过"洞"而不穿过任何三角形的磁通

#### 物理意义

考虑圆柱面，两点 $a$ 和 $b$ 之间的电流定义为：

$$I_{a-b} = \int_a^b (\nabla\chi \times \hat{\mathbf{n}}) \cdot (\hat{\mathbf{n}} \times d\mathbf{l}) = \chi_b - \chi_a$$

如果 $a$ 和 $b$ 位于不同边界，则电势差定义了流过它们之间的电流。对于圆柱面，这对应于环绕圆柱的总电流（角向电流）。

#### 重要规则

**只添加一个洞**：对于圆柱面，虽然有两条边界（顶边和底边），但只需要添加一个"洞"元素。这是因为：
- 角向电流由顶边电势 $\chi_a$ 和底边电势 $\chi_b$ 之差决定
- 如果添加两个洞，会引入**规范冗余**，导致L矩阵奇异

**多连通区域**：对于环面，必须添加两个洞（极向和环向各一个），否则无法支持这两个方向的净电流。

### "闭合"(Closure)元素

#### 目的

对于闭合曲面（如球面或环面），没有边界。由于解只依赖于 $\chi$ 的梯度而不依赖于其绝对值，会引入**规范模糊性**，导致冗余自由度和奇异矩阵。

#### 实现方式

"闭合"元素通过在特定节点固定电势为零来消除规范模糊性：
1. 从活动顶点集合中移除该顶点
2. 该顶点的电势被隐式设为零
3. 系统变为满秩，可以正常求解

#### 代码位置

在 `tw_setup` 子程序中，闭合元素的处理如下：
```fortran
! 将闭合单元转换为顶点并标记
DO i=1,self%nclosures
  ! 选择度数最高的顶点作为闭合点
  ! 从活动顶点映射中移除该顶点
  self%pmap(j) = -i  ! 标记为闭合
END DO
```

### 代码实现

#### 数据结构

```fortran
! 洞元素数据结构
TYPE :: hole_mesh
  INTEGER(i4) :: i0           ! 边界链的起始顶点
  INTEGER(i4) :: n            ! 边界链中的顶点数
  INTEGER(i4), POINTER :: lp(:)  ! 形成闭合边界环的有序顶点列表
  INTEGER(i4), POINTER :: kpc(:) ! CSR格式的单元-顶点链接指针
  INTEGER(i4), POINTER :: lpc(:) ! 每个边界顶点相邻的单元列表
  REAL(r8) :: ptcc(3)         ! 洞边界的几何中心
END TYPE hole_mesh

! 主模型类型
TYPE :: tw_type
  INTEGER(i4) :: nelems       ! 总自由度数 = np_active + nholes + n_vcoils
  INTEGER(i4) :: np_active    ! 活动顶点数（不包括闭合点）
  INTEGER(i4) :: nholes       ! 洞的数量
  INTEGER(i4) :: nclosures    ! 闭合点的数量
  INTEGER(i4), POINTER :: pmap(:)     ! 网格顶点到活动DOF的映射
  INTEGER(i4), POINTER :: kfh(:)      ! 面-洞交互的CSR指针
  INTEGER(i4), POINTER :: lfh(:,:)    ! 面-洞交互列表 (洞索引, 局部边)
  TYPE(hole_mesh), POINTER :: hmesh(:) ! 洞定义数组
END TYPE tw_type
```

#### 溶液向量布局

在时域模拟中，溶液向量 `vec` 的组织如下：

| 索引范围 | 内容 |
|---------|------|
| `1 : np_active` | 网格顶点的标量电势值 |
| `np_active+1 : np_active+nholes` | 洞电流（通过洞的磁通） |
| `np_active+nholes+1 : nelems` | Vcoil（电压驱动线圈）电流 |

#### L矩阵中的洞处理

在计算电感矩阵时，洞自由度的索引计算如下：

```fortran
! 洞DOF索引 = np_active + 洞索引
ik = ABS(row_obj%lfh(1,ii)) + row_obj%np_active

! 洞-洞耦合
jk = ABS(col_obj%lfh(1,jj)) + col_obj%np_active
Lmat(jk,ik) = Lmat(jk,ik) + SIGN(1,row_obj%lfh(1,ii))*SIGN(1,col_obj%lfh(1,jj))*DOT_PRODUCT(...)
```

---

## HODLR矩阵压缩

### 算法原理

电感矩阵 $\mathbf{L}$ 是一个密集矩阵，其元素数量和计算复杂度按 $O(N^2)$ 增长，其中 $N$ 是模型中的元素数量。这限制了大型模型的可扩展性。

**HODLR**（Hierarchical Off-Diagonal Low-Rank，层次对角低秩）方法利用矩阵的结构特性进行近似压缩：

1. **低秩结构**：空间上相距较远的元素之间的相互作用可以用低秩矩阵近似
2. **层次划分**：将网格递归划分为层次结构
3. **压缩存储**：对角块存储为密集矩阵，非对角块存储为低秩近似

这可以将内存和计算复杂度从 $O(N^2)$ 降低到 $O(N \log N)$。

### 空间划分

#### 二叉树划分

ThinCurr使用**二叉树**进行空间划分：

1. 创建包围整个网格的边界框
2. 递归细分：选择位置方差最大的笛卡尔方向进行分割
3. 停止条件：当块中元素数量小于目标大小时停止

#### 块分类

对于每对块 $(i, j)$：

| 分类条件 | 存储方式 |
|---------|---------|
| 对角块 ($i = j$) | 密集矩阵 |
| 近场块（距离近） | SVD压缩 |
| 远场块（距离远） | ACA+压缩 |

分类使用距离-尺寸比：

$$\frac{|\mathbf{c}_i - \mathbf{c}_j|}{r_i + r_j}$$

其中 $\mathbf{c}$ 是块中心，$r$ 是特征尺寸（外接圆半径）。

### 自适应交叉近似(ACA+)

#### 算法概述

**ACA+**（Adaptive Cross Approximation）是一种无需形成完整矩阵即可构建低秩近似的方法：

$$\mathbf{B} \approx \mathbf{U}\mathbf{V}^T = \sum_{k=1}^{r} \mathbf{u}_k \mathbf{v}_k^T$$

#### 算法步骤

1. **初始化**：选择参考行 $i_1$，计算整行元素，找到最大元素位置 $j_1$
2. **迭代**：
   - 选择下一个枢轴行/列（基于残差最大元素）
   - 计算残差行/列（通过核函数即时计算）
   - 更新近似矩阵
3. **收敛判断**：当秩-1更新的范数小于给定容差时停止

#### 关键特性

- **即时计算**：矩阵元素通过积分核即时计算，无需存储完整矩阵
- **自适应选择**：枢轴选择基于残差大小，而非随机采样
- **数值稳定性**：通过参考行列跟踪避免过早终止

#### 收敛准则

算法在以下条件下停止：

$$\|\mathbf{u}_k\|_F \cdot \|\mathbf{v}_k\|_F < \varepsilon$$

其中 $\varepsilon$ 是用户指定的容差。

### 代码实现

#### HODLR操作符类型

```fortran
type, extends(oft_noop_matrix) :: oft_tw_hodlr_op
  ! 层次结构
  INTEGER(4) :: nlevels              ! 层次级别数
  INTEGER(4) :: nblocks              ! 叶子级别的块数
  INTEGER(4) :: ndense               ! 对角（密集）交互数
  INTEGER(4) :: nsparse              ! 压缩的非对角交互数

  ! 容差参数
  REAL(8) :: L_svd_tol               ! L矩阵SVD容差
  REAL(8) :: L_aca_tol               ! L矩阵ACA容差

  ! 存储数组
  TYPE(oft_native_dense_matrix), POINTER :: dense_mats(:)   ! 对角块
  TYPE(oft_native_dense_matrix), POINTER :: aca_U_mats(:)   ! U矩阵
  TYPE(oft_native_dense_matrix), POINTER :: aca_V_mats(:)   ! V矩阵

  ! 层次结构
  TYPE(oft_tw_level), POINTER :: levels(:)  ! 各级别块数组
end type
```

#### 矩阵-向量乘法

```fortran
subroutine tw_hodlr_Lapply(self, a, b)
  ! 1. 应用对角块
  DO i=1,self%ndense
    ! 使用密集矩阵乘法
    b = b + dense_mats(i) * a
  END DO

  ! 2. 应用压缩的非对角块
  DO i=1,self%nsparse
    ! 使用低秩近似: U * V^T * a
    b = b + U_mats(i) * (V_mats(i)^T * a)
  END DO

  ! 3. 应用洞和Vcoil耦合
  b = b + hole_Vcoil_mat * a_holes
end subroutine
```

---

## 时域模拟中的处理

### 时间步进格式

ThinCurr支持两种隐式时间格式：

**后向欧拉**：
$$[\mathbf{L} + \Delta t \cdot \mathbf{R}] \mathbf{I}^{n+1} = \mathbf{L} \cdot \mathbf{I}^n + \Delta t \cdot \mathbf{V}^{n+1}$$

**Crank-Nicolson**：
$$[\mathbf{L} + \frac{\Delta t}{2}\mathbf{R}] \mathbf{I}^{n+1} = [\mathbf{L} - \frac{\Delta t}{2}\mathbf{R}] \mathbf{I}^n + \frac{\Delta t}{2}(\mathbf{V}^{n+1} + \mathbf{V}^n)$$

### 洞电流处理

在时域模拟中，洞电流存储在溶液向量的 `np_active + 1` 到 `np_active + nholes` 位置。这些值代表围绕拓扑非平凡回路流动的净电流。

### 求解器选择

| 模型大小 | 求解方法 |
|---------|---------|
| 小型（<20k元素） | 直接LU分解 |
| 大型（使用HODLR） | 迭代共轭梯度(CG) + 块雅可比预条件 |

### 预条件器

对于HODLR矩阵，使用**块雅可比预条件器**：
- 对角块的LU逆作为预条件
- 保持 $O(N \log N)$ 的求解复杂度

---

## 参考资源

1. **ThinCurr主文档**：`src/docs/ThinCurr/doc_thincurr_main.md`
2. **源代码**：
   - 拓扑处理：`src/physics/thin_wall.F90`
   - HODLR实现：`src/physics/thin_wall_hodlr.F90`
   - 时域求解器：`src/physics/thin_wall_solvers.F90`
3. **示例代码**：
   - Python接口：`examples/ThinCurr/` 目录
   - HODLR示例：`examples/ThinCurr/ports/ports_HODLR.ipynb`
4. **学术论文**：
   - Bebendorf, M. & Rjasanow, S. "Adaptive low-rank approximation of collocation matrices"
   - Erickson, J. & Whittlesey, K. "Greedy optimal homotopy and homology generators"

---

*文档版本：1.0*
*最后更新：2026年3月*