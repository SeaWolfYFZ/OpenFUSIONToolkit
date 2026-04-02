# 66-Hole Mesh Test - Status and Instructions

## Current Status

**问题**: 当前 `tokamak_mesh_EHL2_holes.h5` 有 **65 个 NODESET** (63 个端口 + 2 个同伦基)，跳过了最长的端口边界。

**假设**: 被跳过的端口出现低频振荡电流是因为该边界没有 Hole DOF 稳定。

**测试方案**: 生成包含**所有 64 个端口**作为 Hole 的网格（共 66 个 NODESET），验证是否消除振荡。

## Required File

需要**原始未处理的网格文件** `tokamak_mesh_EHL2.h5`（不包含 NODESET 数据）来重新运行 `ThinCurr_compute_holes.py`。

### 如果已有原始网格：

```bash
# 1. 确保在 example-29 目录
cd /home/yfz/project/OpenFUSIONToolkit-260316/app/ThinCurr-Python-example-29

# 2. 运行修改后的脚本（已修改为包含所有 64 个端口）
python ThinCurr_compute_holes.py \
  --in_file=tokamak_mesh_EHL2.h5 \
  --out_file=tokamak_mesh_EHL2_66holes.h5 \
  --debug

# 3. 验证 NODESET 数量
h5dump -H tokamak_mesh_EHL2_66holes.h5 | grep NUM_NODESETS
# 应显示 66

# 4. 运行测试脚本
python test_66holes.py
```

### 如果没有原始网格：

需要从 CAD 模型重新导出原始网格，或从以下位置查找：
- Tokamak 装置的原始 CAD 文件
- 之前的网格生成脚本输出
- 联系网格生成负责人获取 `tokamak_mesh_EHL2.h5`（无 holes 版本）

## 当前修改

已修改 `ThinCurr_compute_holes.py` 第 659-663 行：

**原代码**（跳过最长边界）:
```python
for k, cycle in enumerate(boundary_cycles[surf_id]):
    if k == cycle_max[1]:
        skipped_holes.append(cycle)  # 跳过最长边界
    else:
        holes.append(cycle)
```

**修改后**（包含所有边界）:
```python
# TEST MODE: Include ALL 64 port boundaries as holes
for cycle in boundary_cycles[surf_id]:
    holes.append(cycle)  # Include ALL boundaries
```

## 预期结果

### 如果 66-Hole 方案有效：
1. 所有 64 个端口电流分布对称
2. 被跳过端口的振荡消失
3. L 矩阵条件数在可接受范围 (< 10¹²)

### 如果 66-Hole 方案失败：
1. L 矩阵奇异或条件数过大 (> 10¹⁴)
2. 求解器不收敛
3. 需要考虑其他方案（如添加边界稳定项）

## 数学分析

### 拓扑正确性
- **理论独立 Hole 数**: 65 (63 端口边界 + 2 同伦基)
- **实际创建 Hole 数**: 66 (64 端口边界 + 2 同伦基)
- **冗余度**: 1

### 数值可行性
虽然拓扑上有 1 个冗余，但由于：
1. 64 个端口在空间上分离
2. 互感矩阵 L_hh 的对角项提供数值稳定性
3. 几何分离使矩阵保持良态

**预计条件数**: 10⁸ - 10¹⁰（可接受范围）

## 下一步

1. **获取原始网格** `tokamak_mesh_EHL2.h5`
2. **运行测试**生成 66-hole 网格
3. **比较结果**与 65-hole 情况
4. **决定方案**: 66-hole 或其他稳定化方法

## 联系信息

如有问题，请检查：
- `test_66holes.py` - 测试脚本
- `ThinCurr_compute_holes.py` (第 659-663 行) - 修改位置
- `reference/ThinCurr-doc.md` - 拓扑处理详细文档
