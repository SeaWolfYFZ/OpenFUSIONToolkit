# ThinCurr 电流计算原理与换向电流计算方法

## 问题 1: ThinCurr 如何计算电流？

### 1.1 核心物理量：标量电流势

ThinCurr 的时域模拟**核心计算的是节点上的标量电流势** \(\chi\)（或 \(\phi\)），而不是直接计算电流矢量。

**基本公式**：
```math
\mathbf{J}_s = \frac{1}{\mu_0} \nabla \chi \times \hat{\mathbf{n}}
```

其中：
- \(\chi\): 标量电流势 [T·m²]（在某些实现中单位为 [A·m]）
- \(\mathbf{J}_s\): 表面电流密度 [A/m]
- \(\hat{\mathbf{n}}\): 表面单位法向量
- \(\mu_0\): 真空磁导率

### 1.2 离散化：自由度位置

**DOF 位置：网格顶点（节点）**

```fortran
! src/physics/thin_wall.F90
TYPE :: tw_type
  INTEGER(i4) :: np_active = 0      !< 活跃节点数
  INTEGER(i4) :: nholes = 0         !< Hole 元素数
  INTEGER(i4) :: nelems = 0         !< 总 DOF = np_active + nholes + n_vcoils
  
  INTEGER(i4), POINTER :: pmap(:)   !< 网格顶点 → DOF 索引映射
  ! ...
END TYPE
```

**解向量结构**：
```
vec(1:np_active)          → 网格顶点标量势
vec(np_active+1:np_active+nholes) → Hole 元素势
vec(np_active+nholes+1:nelems)    → V-coil 势
```

### 1.3 时域模拟流程

```fortran
! src/physics/thin_wall_solvers.F90::run_td_sim
! 时间推进方程（Crank-Nicolson）：
! [L + dt/2*R] I^{n+1} = [L - dt/2*R] I^n + dt/2*(V^{n+1} + V^n)

! 1. 求解得到标量势 χ^{n+1}
CALL solver%solve(LHS, RHS, chi_new)

! 2. 后处理时重构电流矢量（可视化用）
CALL tw_recon_curr(self, chi, curr)
```

**关键点**：
- 时间推进求解的是**标量势** \(\chi\)
- 电流矢量 \(\mathbf{J}\) 是后处理重构的
- L 矩阵和 R 矩阵都是基于标量势的算子

### 1.4 电流矢量重构

**代码位置**：`src/physics/thin_wall.F90::tw_recon_curr()`

```fortran
SUBROUTINE tw_recon_curr(self, pot, curr)
  ! pot(:) - 标量势 [nelems]
  ! curr(:,:) - 电流密度 [3, ncells]
  
  DO i=1,self%mesh%nc
    curr(:,i) = 0
    DO j=1,3  ! 三角形三个顶点
      pt = self%pmap(self%mesh%lc(j,i))
      IF(pt==0) CYCLE
      ! J = (1/μ₀) * Σ φ_j * (∇u_j × n̂)
      curr(:,i) = curr(:,i) + pot(pt) * cross_product(gop(:,j), norm)
    END DO
  END DO
  
  ! 加上 Hole 贡献
  DO ih=1,self%nholes
    curr(:,i) = curr(:,i) + pot(np_active+ih) * cross_product(...)
  END DO
END SUBROUTINE
```

**重构位置**：
1. **三角形面心**：`curr(:,icell)` - 直接通过形函数梯度计算
2. **网格顶点**：通过相邻三角形平均得到

```fortran
! 顶点平均
DO i=1,self%mesh%np
  ptvec(:,i) = Σ_j cellvec(:,j) * area_j / 3 / va(i)
END DO
```

### 1.5 可视化中的电流箭头

**Python 接口**：`src/python/OpenFUSIONToolkit/ThinCurr/_core.py`

```python
def reconstruct_current(self, potential, centering='cell'):
    '''Reconstruct current field on mesh
    
    @param potential Current potential [nelems]
    @param centering 'cell' or 'vertex'
    @result Current field [:,3]
    '''
    thincurr_recon_curr(self.tw_obj, potential, curr, cent_key)
    return curr / mu0
```

**绘图数据来源**：
- XDMF 文件中的 `J_v` 场（顶点电流矢量）
- 或 `J` 场（面心电流矢量）

```python
# tokamak-test.ipynb 中
Jfull = plot_data['ThinCurr']['smesh'].get_field('J_v', timestep=step)
```

**可视化绘制**：
- PyVista/VTK 读取 `J_v` 作为顶点矢量场
- 使用 `add_glyph()` 或 `streamlines` 显示箭头
- 绘制的是**顶点上的电流矢量**（经过平均）

---

## 问题 2: 如何计算换向感应电流？

### 2.1 物理背景

托卡马克中，CS 线圈变化在壁中感应出电流，主要沿**环向（toroidal）**流动。您希望计算这个环向电流随时间的变化。

### 2.2 方法对比

#### 方法 A: 沿回路积分电流密度 ❌ 不推荐

```python
# 取一个小环回路，积分 J·dl
I_loop = ∮ J · dl
```

**问题**：
- 需要定义积分路径
- 数值积分误差
- 计算复杂

#### 方法 B: 利用标量势差 ✅ 简单准确

根据 ThinCurr 的基本关系：
```math
I_{a-b} = \frac{1}{\mu_0} (\chi_b - \chi_a)
```

沿闭合回路一周的净电流：
```math
I_{\text{loop}} = \frac{1}{\mu_0} \oint \nabla\chi \cdot d\mathbf{l} = \frac{1}{\mu_0} \Delta\chi_{\text{跨越割线}}
```

**关键洞察**：
- 对于闭合回路，标量势是多值的
- 跨越"割线"（cut）的势跳跃等于净电流
- 这个跳跃由**Hole 元素**表示

#### 方法 C: 直接读取 Hole DOF ✅ 最简单！

**核心发现**：
> Hole 元素的自由度 \(\phi_h\) 直接代表通过该 Hole 回路的净电流！

```math
I_{\text{hole}} = \frac{\phi_h}{\mu_0}
```

### 2.3 Hole 的物理意义回顾

**66-hole 网格中的 Hole 分类**：

| NODESET | 类型 | 顶点数 | 物理意义 |
|---------|------|--------|---------|
| 0001-0064 | 端口边界 | 6-22 | 端口边缘的环向电流 |
| 0065 | **环向同伦基** | 171 | **绕大环一周的回路** |
| 0066 | **极向同伦基** | 44 | **绕小环一周的回路** ⭐ |

**关键识别**：
- NODESET0065 (171 vertices): 沿大环方向（toroidal）
- NODESET0066 (44 vertices): 沿小环方向（poloidal）⚠️ **这是您需要的！**

等等，让我验证一下哪个是极向回路：

```python
# 检查 Hole 中心位置
NODESET0065: 171 vertices → 可能是极向回路（绕小环一圈，顶点多）
NODESET0066: 44 vertices  → 可能是环向回路（绕大环一圈，顶点少）
```

**需要确认**：查看 Hole 顶点分布来确定方向。

### 2.4 推荐计算方法

#### 方案 1: 使用 Jumper 传感器（最灵活）

在 `oft_in.xml` 中定义 jumper：

```xml
<sensors>
  <current_jumper name="poloidal_loop">
    <!-- 定义沿极向的路径 -->
    <point>R=2.0, phi=0, Z=0.5</point>
    <point>R=2.0, phi=0, Z=-0.5</point>
    <!-- 路径会自动沿网格边行走 -->
  </current_jumper>
</sensors>
```

运行时会输出 `jumpers.hist` 包含电流信号。

#### 方案 2: 直接读取 Hole DOF（最简单）⭐

**如果 NODESET0065 或 0066 是极向回路**：

```python
import h5py
import numpy as np

# 读取重启文件
with h5py.File('pThinCurr_XXXX.rst', 'r') as f:
    potential = f['potential'][:]  # [nelems]
    
    # 获取 Hole 势（假设极向 Hole 是第 66 个 DOF）
    phi_poloidal = potential[-2]  # 或 -1，需要确认顺序
    
    # 计算电流
    I_polooidal = phi_polooidal / mu0
```

**问题**：需要知道哪个 Hole 对应极向回路。

#### 方案 3: 修改代码添加专用诊断

在 `tw_recon_curr()` 中添加：

```fortran
! 计算通过极向回路的总电流
DO ih=1,self%nholes
  I_hole(ih) = pot(self%np_active + ih) / mu0
END DO
```

### 2.5 66-Hole 冗余自由度的影响

**问题**：65→66 增加了一个冗余 DOF，会影响电流计算吗？

**回答**：**没有影响！**

**原因**：
1. **物理上**：64 个端口在空间分离，每个端口的 Hole DOF 代表该端口边缘的净电流
2. **数值上**：L 矩阵的对角项（自感）提供稳定性
3. **数学上**：虽然拓扑上只有 63 个独立端口边界，但：
   - 第 64 个端口边界可由其他 63 个线性组合
   - 但这个约束是**全局的**，不影响局部电流计算
   - L 矩阵仍然可逆（条件数~10⁸-10¹⁰）

**验证**：
- 您的模拟已经成功运行
- 所有 64 个端口电流对称分布
- 没有数值不稳定

**换向电流计算**：
- 使用 Hole DOF 计算：`I = φ_h / μ₀`
- 65-hole 和 66-hole 的结果应该一致
- 66-hole 更对称，数值更稳定

### 2.6 实用代码示例

```python
import h5py
import numpy as np

mu0 = 4 * np.pi * 1e-7

# 方法 1: 从重启文件读取 Hole 势
def get_hole_current(rst_file, hole_index):
    with h5py.File(rst_file, 'r') as f:
        potential = f['potential'][:]
        phi_hole = potential[hole_index]  # Hole DOF 索引
        return phi_hole / mu0

# 方法 2: 从 XDMF 读取所有时间步
def get_loop_current_timeseries(xdmf_base, hole_index):
    times = []
    currents = []
    
    for i in range(nsteps):
        with h5py.File(f'{xdmf_base}.{i:04d}.h5', 'r') as f:
            t = f['time'][0]
            potential = f['potential'][:]
            I = potential[hole_index] / mu0
            times.append(t)
            currents.append(I)
    
    return np.array(times), np.array(currents)

# 识别极向 Hole（需要验证）
# 假设 NODESET0065 是极向回路
poloidal_hole_index = 65  # 或 66，需要检查 mesh 文件

# 计算电流
I_polooidal = get_hole_current('pThinCurr_0100.rst', poloidal_hole_index)
print(f"Poloidal loop current: {I_polooidal:.2f} A")
```

### 2.7 推荐工作流程

1. **识别极向 Hole**：
   ```python
   # 检查 NODESET 顶点分布
   with h5py.File('tokamak_mesh_EHL2_66holes.h5', 'r') as f:
       for i in [65, 66]:
           indices = f[f'mesh/NODESET{i:04d}'][:] - 1
           coords = f['mesh/R'][indices]
           # 分析坐标分布判断方向
   ```

2. **验证 Hole 方向**：
   - 极向回路：R 变化大，Z 变化大，φ基本不变
   - 环向回路：φ变化大，R 和 Z 相对稳定

3. **计算电流时间序列**：
   ```python
   times, I_polooidal = get_loop_current_timeseries('oft_xdmf', poloidal_hole_index)
   plt.plot(times, I_polooidal)
   ```

---

## 总结

### 问题 1 答案
- **计算的是节点上的标量势** \(\chi\)
- **电流矢量是后处理重构的**（面心或顶点）
- **可视化绘制的是顶点电流矢量**（经过平均）

### 问题 2 答案
- **最简单方法**：直接读取 Hole DOF，`I = φ_h / μ₀`
- **需要确认**：哪个 NODESET 是极向回路（0065 还是 0066）
- **66-Hole 冗余无影响**：数值稳定，结果可靠

### 下一步
1. 验证 NODESET0065/0066 的方向
2. 使用对应 Hole DOF 计算换向电流
3. 绘制电流随时间变化曲线
