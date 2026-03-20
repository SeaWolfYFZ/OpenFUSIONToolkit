# ThinCurr 网格拓扑与 HODLR 实现详解

> 本文档面向等离子体物理学专业的学生和研究人员，详细解释 ThinCurr 模块中的网格拓扑处理（Hole 和 Closure 元素）以及 HODLR 矩阵压缩算法的实现原理。

[TOC]

## 第一部分：网格拓扑处理

### 1.1 为什么需要拓扑处理？

ThinCurr 使用**标量电流势**（scalar current potential）\(\chi\) 来表示薄壁结构中的表面电流：

\[
\mathbf{J}_s = \nabla \chi \times \hat{\mathbf{n}}
\]

其中 \(\hat{\mathbf{n}}\) 是表面的单位法向量。这种表示方法自动满足电流的无散条件（\(\nabla \cdot \mathbf{J}_s = 0\)），但在处理**多连通区域**（multiply-connected domains）时会遇到问题。

#### 1.1.1 拓扑学基础概念

在拓扑学中，曲面的**连通性**由其"洞"的数量（即贝蒂数/Betti number）描述：

| 几何形状 | 拓扑类型 | 贝蒂数 (b₁) | 所需 Hole 数 |
|---------|---------|------------|-------------|
| 圆盘（disk） | 单连通 | 0 | 0 |
| 圆柱（cylinder） | 多连通 | 1 | 1 |
| 环面（torus） | 多连通 | 2 | 2 |
| 球面（sphere） | 闭曲面 | 0 | 0（但需要 Closure） |

**贝蒂数 b₁** 表示曲面上线性无关的闭合回路的数量。这些回路无法通过连续变形收缩为一个点。

### 1.2 Hole 元素：处理多连通几何

#### 1.2.1 Hole 的物理意义

考虑圆柱面上的电流。如果只用单值标量势 \(\chi\)，则沿圆周方向积分：

\[
I_{\text{环向}} = \oint \nabla \chi \cdot d\mathbf{l} = \chi_{\text{终点}} - \chi_{\text{起点}} = 0
\]

因为起点和终点重合，单值势的差为零！这显然无法描述环向电流。

**解决方案**：引入 Hole 元素，允许势在跨越 Hole 定义的"割线"（cut）时发生跳跃，从而使 \(\chi\) 成为**多值函数**。

数学上，跨越割线的势跳跃量等于通过该回路的净电流：

\[
\Delta\chi_{\text{hole}} = \chi_+ - \chi_- = \mu_0 I_{\text{通过回路}}
\]

#### 1.2.2 Hole 的数学实现

在 ThinCurr 中，每个 Hole 对应：
1. **一个闭合的边界顶点链**：定义 Hole 的几何位置
2. **一个额外的自由度**：Hole 上的常数势值 \(\phi_h\)
3. **面 -Hole 耦合**：相邻三角形与 Hole 的相互作用

总自由度数：
\[
N_{\text{DOF}} = N_{\text{active}} + N_{\text{holes}} + N_{\text{vcoils}}
\]

其中：
- \(N_{\text{active}}\)：活跃的内部/边界节点（排除 closure）
- \(N_{\text{holes}}\)：Hole 元素数量
- \(N_{\text{vcoils}}\)：电压驱动的线圈数量

#### 1.2.3 Hole 的自动生成算法

ThinCurr 使用 **Greedy Homotopy 算法**（基于 Erickson & Whittlesey 的方法）自动生成 Hole：

```python
# ThinCurr_compute_holes.py 中的核心流程

1. 边界识别
   - 识别所有边界循环（boundary cycles）
   - 除最长边界外，其余都加入 Hole 列表

2. 曲面封闭
   - 对每个边界循环，添加"扇形"三角形进行人工封闭
   - 封闭边赋予高权重，避免拓扑回路穿过

3. 贪婪同伦基搜索
   - 从基点（basepoint）出发，使用 Dijkstra 算法构建最短路径树
   - 在对偶图中找出生成树
   - 剩余边对应同伦基的生成元

4. 优化
   - 从不同基点重新运行搜索
   - 选择最短的回路组合作为最终 Hole 定义
```

**关键代码位置**：
- `src/utilities/scripts/ThinCurr_compute_holes.py`: 主脚本
- `compute_greedy_homotopy_basis()`: 核心算法函数
- `fixup_loop()`: 回路修复和优化

#### 1.2.4 Hole 在矩阵组装中的处理

在电感矩阵 L 的组装中，Hole 元素通过**面-Hole 耦合**与普通三角形 DOF 相互作用：

```fortran
! src/physics/thin_wall.F90 中的 face-hole 连通性
INTEGER(i4), POINTER :: kfh(:)  ! CSR 格式指针数组 [mesh%nc+1]
INTEGER(i4), POINTER :: lfh(:,:) ! 面 -Hole 连接列表 [2, nfh]
! lfh(1,ii) = 带符号的 Hole 索引
! lfh(2,ii) = 三角形上的局部边索引 (1,2,3)
```

电感矩阵块的结构：
```
L = [ L_aa   L_ah ]
    [ L_ha   L_hh ]

其中：
- L_aa: 活跃节点之间的互感 [np_active × np_active]
- L_ah: 活跃节点与 Hole 的耦合 [np_active × nholes]
- L_hh: Hole 之间的自感/互感 [nholes × nholes]
```

### 1.3 Closure 元素：规范固定

#### 1.3.1 为什么需要 Closure？

对于**闭曲面**（如完整环面、球面），没有自然边界来固定电势参考点。此时电感矩阵 L 存在**规范自由度**：

\[
\chi \rightarrow \chi + C \quad \Rightarrow \quad \mathbf{J}_s \text{不变}
\]

这导致 L 矩阵奇异（不满秩），无法求逆。

#### 1.3.2 Closure 的实现

Closure 通过**固定某点的势为零**来消除规范自由度：

1. 选择一个三角形（通常是网格中的某个单元）
2. 从该三角形的三个顶点中，选择连接数最多的顶点
3. 将该顶点的 DOF 从系统中移除（pmap = 0）

```fortran
! src/physics/thin_wall.F90 中的 Closure 处理
DO i=1,self%nclosures
  ! 选择连接最多的顶点作为 closure 点
  DO k=1,3
    IF(self%mesh%kpc(self%mesh%lc(k,self%closures(i))+1) - &
       self%mesh%kpc(self%mesh%lc(k,self%closures(i))) > l) THEN
      l = ...
      j = k
    END IF
  END DO
  j = self%mesh%lc(j,self%closures(i))
  self%pmap(j) = -i  ! 临时标记为 closure
END DO
```

#### 1.3.3 Hole 与 Closure 的区别

| 特性 | Hole | Closure |
|-----|------|---------|
| **目的** | 支持多连通几何的净电流 | 固定闭曲面的电势规范 |
| **数量** | 等于贝蒂数 b₁ | 等于包围的体积数（通常为 1） |
| **DOF** | 增加自由度 | 移除自由度 |
| **物理意义** | 拓扑非平凡回路 | 规范固定点 |
| **矩阵影响** | 增加 L 矩阵维度 | 确保 L 矩阵满秩 |

### 1.4 时域模拟中的拓扑处理

在时域模拟（`run_td_sim`）中，Hole 和 Closure 的处理体现在：

#### 1.4.1 时间推进方程

使用 Crank-Nicolson 或向后 Euler 格式：

\[
\left[\mathbf{L} + \frac{\Delta t}{2}\mathbf{R}\right] \mathbf{I}^{n+1} = 
\left[\mathbf{L} - \frac{\Delta t}{2}\mathbf{R}\right] \mathbf{I}^{n} + 
\frac{\mathbf{V}^{n+1} + \mathbf{V}^{n}}{2}
\]

其中电流向量包含所有 DOF：
\[
\mathbf{I} = [\mathbf{I}_{\text{active}}, \mathbf{I}_{\text{holes}}, \mathbf{I}_{\text{vcoils}}]^T
\]

#### 1.4.2 电流跳变传感器（Jumper Sensors）

Jumper 传感器测量沿路径的净电流，需要包含 Hole 的贡献：

```fortran
! src/physics/thin_wall_solvers.F90 中的 jumper 计算
tmp = 0.d0
val_prev = 0.d0
! 沿路径积分势差
DO k=1,sensors%jumpers(j)%np-1
  ind1 = self%pmap(sensors%jumpers(j)%points(k+1))
  IF(ind1 > 0) THEN
    tmp = tmp + vals(ind1) - val_prev
    val_prev = vals(ind1)
  END IF
END DO
! 加上 Hole 的贡献
DO k=1,self%nholes
  tmp = tmp + vals(self%np_active+k) * sensors%jumpers(j)%hole_facs(k)
END DO
jumpout(j+1) = tmp / mu0
```

其中 `hole_facs(k)` 是第 k 个 Hole 对该 Jumper 的耦合系数，由拓扑决定。

---

## 第二部分：HODLR 矩阵压缩算法

### 2.1 为什么需要 HODLR？

传统边界元方法（BEM）产生的电感矩阵 L 是**稠密矩阵**，存储和计算复杂度均为 \(O(N^2)\)：

| 网格规模 | 内存需求 | 计算时间 |
|---------|---------|---------|
| 1,000 DOF | ~8 MB | ~1 秒 |
| 10,000 DOF | ~800 MB | ~100 秒 |
| 100,000 DOF | ~80 GB | ~3 小时 |
| 1,000,000 DOF | ~8 TB | 不可行 |

对于实际聚变装置（如 ITER、CFETR），网格规模往往超过 10 万 DOF，必须使用压缩算法。

### 2.2 HODLR 的基本思想

#### 2.2.1 低秩结构

电感矩阵的**远场相互作用**具有低秩特性：物理上分离较远的区域之间，磁场耦合可以用少数几个"模式"近似描述。

数学上，对于空间分离的两个区域 i 和 j，对应的矩阵块 \(L_{ij}\) 可以近似为：

\[
L_{ij} \approx U \cdot V^T
\]

其中 \(U \in \mathbb{R}^{m \times k}\)，\(V \in \mathbb{R}^{n \times k}\)，且 \(k \ll \min(m,n)\)。

#### 2.2.2 层次化分块

HODLR 的核心是**递归二分**的空间划分：

```
Level 1: [              全区域              ]  (1 块)
Level 2: [    左半部    ][    右半部    ]  (2 块)
Level 3: [ 左 ][ 左 ][ 右 ][ 右 ]  (4 块)
...
Level L: [叶节点块，每块~1500 DOF]  (N/1500 块)
```

每个层级上，矩阵分块分为：
- **对角块**（Near-field）：相邻区域，稠密存储
- **非对角块**（Far-field）：远离区域，低秩压缩

### 2.3 HODLR 实现流程

#### 2.3.1 空间划分（`tw_hodlr_setup`）

```fortran
! src/physics/thin_wall_hodlr.F90 中的划分逻辑

SUBROUTINE tw_hodlr_setup(self, required)
  1. 读取配置参数（thincurr_hodlr_options）
     - target_size = 1500  (叶节点目标大小)
     - L_svd_tol = 1e-6    (SVD 截断容差)
     - L_aca_rel_tol = 0.1 (ACA 相对容差)
  
  2. 递归二分网格
     - 计算当前块的包围盒
     - 沿标准差最大的坐标轴切割
     - 直到块大小 <= target_size
  
  3. 分类块间相互作用
     - 计算块中心距离 / (块半径之和)
     - 比值 > 1.1 → 远场 (ACA 压缩)
     - 比值 ≤ 1.1 → 近场 (稠密积分)
  
  4. 建立掩码矩阵
     - 避免在不同层级重复计算
     - mat_mask(j,k) = -1 表示已被父层级覆盖
END SUBROUTINE
```

#### 2.3.2 矩阵组装（`tw_hodlr_Lcompute`）

```fortran
SUBROUTINE tw_hodlr_Lcompute(self, save_file)
  1. Hole/V-coil 耦合矩阵（稠密）
     CALL tw_compute_LmatHole(...)
  
  2. 对角块（稠密积分）
     DO i = 1, ndense
       CALL tw_compute_Lmatblock(...)  ! 直接 Biot-Savart 积分
     END DO
  
  3. 非对角块（ACA+ 压缩）
     DO i = 1, nsparse
       IF (远场) THEN
         CALL aca_approx(...)  ! ACA 迭代
         CALL compress_block(...)  ! SVD 后压缩
       ELSE
         CALL tw_compute_Lmatblock(...)  ! 回退到稠密
       END IF
     END DO
END SUBROUTINE
```

### 2.4 ACA+ 算法详解

#### 2.4.1 算法原理

Adaptive Cross Approximation (ACA+) 通过**部分 pivoting**迭代构建低秩近似：

**目标**：近似矩阵块 \(A \in \mathbb{R}^{m \times n}\) 为 \(A \approx U_k V_k^T\)，其中 \(U_k \in \mathbb{R}^{m \times k}\)，\(V_k \in \mathbb{R}^{k \times n}\)。

**迭代步骤**：

```
初始化：k = 1, R_0 = A, 选择随机行索引 i_1

迭代 k = 1, 2, ...:
  1. 选择 pivot 行：i_k = argmax |(R_{k-1})_{i,:}|
  2. 计算残差行：v_k = R_{k-1}(i_k, :)
  3. 选择 pivot 列：j_k = argmax |(v_k)_j|
  4. 归一化：v_k = v_k / v_k(j_k)
  5. 计算残差列：u_k = R_{k-1}(:, j_k)
  6. 更新残差：R_k = R_k - u_k * v_k^T
  7. 如果 ||u_k|| * ||v_k|| < ε * ||A_k||，停止
  8. k = k + 1
```

**关键优势**：只需计算矩阵的**部分行和列**，无需组装完整矩阵块。

#### 2.4.2 ACA+ 的 ThinCurr 实现

```fortran
! src/physics/thin_wall_hodlr.F90 中的 ACA+ 实现

SUBROUTINE aca_approx(isparse, tol, size_out)
  ! 初始化
  ALLOCATE(US(N, max_iter), VS(max_iter, M))
  Iref = 随机选择参考行
  Jref = 随机选择参考列
  
  ! 计算初始参考行/列（动态调用 Biot-Savart 积分）
  CALL tw_compute_Lmatblock(..., RIref, ...)  ! 第 Iref 行
  CALL tw_compute_Lmatblock(..., RJref, ...)  ! 第 Jref 列
  
  DO k = 1, max_iter
    ! 1. 选择 pivot 行（排除已选行）
    Jstar = max_masked(ABS(RIref), prevJstar)
    
    ! 2. 选择 pivot 列（排除已选列）
    Istar = max_masked(ABS(RJref), prevIstar)
    
    ! 3. 决定先 pivoting 行还是列
    IF (|RIref(Jstar)| > |RJref(Istar)|) THEN
      ! 先行 pivoting
      CALL compute_row(Jstar, ...)
      CALL compute_col(Jstar, ...)
    ELSE
      ! 先列 pivoting
      CALL compute_col(Istar, ...)
      CALL compute_row(Istar, ...)
    END IF
    
    ! 4. 更新残差
    RIref = RIref - u_k(Iref) * v_k
    RJref = RJref - u_k * v_k(Jref)
    
    ! 5. 检查收敛（使用快速 Frobenius 范数递推）
    norm_Sk = SQRT(norm_Sk_prev + 2*dot_product(...) + norm_u*norm_v)
    IF (norm_u * norm_v < tol * norm_Sk) EXIT
  END DO
  
  ! 6. SVD 后压缩
  CALL compress_block(...)
END SUBROUTINE
```

#### 2.4.3 快速 Frobenius 范数计算

直接计算 \(\|S_k\|_F\) 需要 \(O(mn)\) 操作，破坏 ACA 的效率。ThinCurr 使用**递推公式**：

\[
\|S_k\|_F^2 = \|S_{k-1}\|_F^2 + 2\sum_{j=1}^{k-1} (u_k^T u_j)(v_j^T v_k) + \|u_k\|_F^2 \|v_k\|_F^2
\]

复杂度：\(O(k(m+n))\)，保持 ACA 的线性复杂度。

### 2.5 HODLR 在时域模拟中的应用

#### 2.5.1 矩阵 - 向量乘法

在时间推进中，需要计算 \(\mathbf{L} \cdot \mathbf{I}\)：

```fortran
! src/physics/thin_wall_hodlr.F90 中的 HODLR 矩阵乘法

SUBROUTINE tw_hodlr_Lapply(self, x, y)
  ! 1. 对角块贡献（稠密矩阵乘法）
  DO i = 1, ndense
    level = dense_blocks(1,i)
    j = dense_blocks(2,i)
    k = dense_blocks(3,i)
    y(block_j) = y(block_j) + MATMUL(dense_mats(i), x(block_k))
  END DO
  
  ! 2. 非对角块贡献（低秩矩阵乘法）
  DO i = 1, nsparse
    level = sparse_blocks(1,i)
    j = sparse_blocks(2,i)
    k = sparse_blocks(3,i)
    ! y += U * (V^T * x)  两次矩阵 - 向量乘法
    tmp = MATMUL(TRANSPOSE(aca_V_mats(i)), x(block_k))
    y(block_j) = y(block_j) + MATMUL(aca_U_mats(i), tmp)
  END DO
  
  ! 3. Hole/V-coil 贡献
  IF (nholes > 0) THEN
    y(active) = y(active) + MATMUL(hole_Vcoil_mat, x(holes))
    y(holes) = y(holes) + MATMUL(TRANSPOSE(hole_Vcoil_mat), x(active))
  END IF
END SUBROUTINE
```

**复杂度分析**：
- 对角块：\(O(N_{\text{leaf}}^2 \times N_{\text{blocks}}) = O(N \times N_{\text{leaf}})\)
- 非对角块：\(O(k \times N)\)，其中 k 为平均秩
- 总体：\(O(N \log N)\)

#### 2.5.2 预处理技术

对于大规模系统，使用**块-Jacobi 预处理**加速迭代求解：

```fortran
! src/physics/thin_wall_hodlr.F90 中的块-Jacobi 预处理

TYPE :: oft_tw_hodlr_rbjpre
  TYPE(oft_tw_hodlr_op), POINTER :: mf_obj  ! HODLR 矩阵
  TYPE(oft_native_matrix), POINTER :: Rmat  ! 电阻矩阵
  TYPE(rmat_container), POINTER :: inverse_mats  ! 对角块逆矩阵
  
CONTAINS
  PROCEDURE :: apply => rbjprecond_apply
END TYPE

SUBROUTINE rbjprecond_apply(self, x, y)
  ! 对每个对角块求解局部系统
  DO i = 1, nblocks
    ! (α*L_diag + β*R_diag) * y_i = x_i
    y(block_i) = MATMUL(inverse_mats(i), x(block_i))
  END DO
END SUBROUTINE
```

预处理矩阵由对角块的 LU 分解预先计算，每次迭代只需前代/回代。

### 2.6 性能与精度权衡

#### 2.6.1 压缩参数选择

| 参数 | 推荐值 | 影响 |
|-----|-------|------|
| `target_size` | 1000-2000 | 较小→树更深，压缩更好但 overhead 更高 |
| `L_svd_tol` | 1e-6 - 1e-8 | 较小→精度更高，秩更大 |
| `L_aca_rel_tol` | 0.01 - 0.1 | 较小→ACA 更精确，迭代更多 |
| `aca_min_its` | 20-30 | 确保 ACA 收敛稳定性 |

#### 2.6.2 压缩效果示例

对于 tokamak 真空室网格（~40,000 DOF）：

```
Partitioning grid for block low rank compressed operators
  nBlocks =                  32
  Avg block size =          618
  # of SVD =                210
  # of ACA =                166

Building block low rank inductance operator
  Compression ratio:   8.3%  (3.23E+07 / 3.91E+08)
  Time = 24s
```

**内存节省**：12 倍
**计算加速**：~10 倍（矩阵组装）+ ~5 倍（矩阵 - 向量乘法）

---

## 第三部分：实践指南

### 3.1 生成 Hole 和 Closure

使用提供的 Python 脚本自动生成拓扑元素：

```bash
python ThinCurr_compute_holes.py \
  --in_file=tokamak_mesh.h5 \
  --out_file=tokamak_mesh_holes.h5 \
  --plot_final \
  --optimize_holes
```

**常用选项**：
- `--plot_final`: 可视化最终的同伦基
- `--optimize_holes`: 尝试优化 Hole 路径（更短、更平滑）
- `--ref_point x y z`: 指定同伦基的基点位置
- `--debug`: 输出详细调试信息

### 3.2 配置 HODLR 参数

在 XML 输入文件中设置：

```xml
<oft>
  <thincurr>
    <!-- 物理参数 -->
    <eta>3.7e-5, 2.467e-5</eta>
    
    <!-- HODLR 参数 -->
    <thincurr_hodlr_options>
      <target_size>1500</target_size>
      <aca_min_its>20</aca_min_its>
      <L_svd_tol>1.0e-6</L_svd_tol>
      <L_aca_rel_tol>0.05</L_aca_rel_tol>
      <B_svd_tol>1.0e-2</B_svd_tol>
      <B_aca_rel_tol>0.1</B_aca_rel_tol>
    </thincurr_hodlr_options>
  </thincurr>
</oft>
```

### 3.3 Python 接口示例

```python
from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr

# 初始化环境
myOFT = OFT_env(nthreads=28)
tw = ThinCurr(myOFT)

# 设置模型（包含自动加载 holes）
tw.setup_model(
    mesh_file='tokamak_mesh_holes.h5',
    xml_filename='oft_in.xml'
)

# 计算 HODLR 压缩的电感矩阵
tw.compute_Lmat(
    use_hodlr=True,
    cache_file='L_matrix.h5'  # 缓存以加速后续运行
)

# 时域模拟
dt = 1.0e-3  # 时间步长 [s]
nsteps = 1000  # 步数
tw.run_td(dt, nsteps, status_freq=10)
```

### 3.4 调试技巧

1. **检查 Hole 定义**：
   ```python
   print(f"Number of holes: {tw.nholes}")
   print(f"Number of closures: {tw.nclosures}")
   ```

2. **验证 HODLR 精度**：
   ```python
   # 比较 HODLR 与稠密矩阵的矩阵 - 向量乘积
   y_hodlr = hodlr_L @ x
   y_dense = dense_L @ x
   error = np.linalg.norm(y_hodlr - y_dense) / np.linalg.norm(y_dense)
   print(f"Relative error: {error:.2e}")
   ```

3. **监控 ACA 收敛**：
   在 `thincurr_hodlr_options` 中设置较小的 `aca_min_its` 并观察警告信息。

---

## 参考文献

1. **拓扑学基础**:
   - Erickson, J., & Whittlesey, K. (2005). "Greedy optimal homotopy and homology generators." SODA.
   - Hatcher, A. (2002). "Algebraic Topology." Cambridge University Press.

2. **HODLR 与 ACA**:
   - Bebendorf, M., & Rjasanow, S. (2003). "Adaptive low-rank approximation of collocation matrices." Computing.
   - Börm, S. (2010). "Efficient Numerical Methods for Non-local Operators." EMS Press.

3. **ThinCurr 文档**:
   - OpenFUSIONToolkit 官方文档：`src/docs/ThinCurr/doc_thincurr_main.md`
   - 代码注释：`src/physics/thin_wall*.F90`

---

## 附录：关键代码位置索引

| 功能 | 文件 | 函数/子程序 |
|-----|------|------------|
| Hole 拓扑分析 | `ThinCurr_compute_holes.py` | `compute_greedy_homotopy_basis()` |
| Hole 设置 | `thin_wall.F90` | `tw_setup()`, `get_hole_pseq()` |
| HODLR 设置 | `thin_wall_hodlr.F90` | `tw_hodlr_setup()` |
| HODLR 矩阵组装 | `thin_wall_hodlr.F90` | `tw_hodlr_Lcompute()` |
| ACA+ 算法 | `thin_wall_hodlr.F90` | `aca_approx()` |
| SVD 压缩 | `thin_wall_hodlr.F90` | `compress_block()` |
| 时域模拟 | `thin_wall_solvers.F90` | `run_td_sim()` |
| Jumper 计算 | `thin_wall_solvers.F90` | `run_td_sim()` 内部 |
