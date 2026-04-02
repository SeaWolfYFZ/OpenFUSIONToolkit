# ==============================================================================
# 7. pyvista可视化检验
# ==============================================================================

import pyvista as pv
import numpy as np
import h5py

hdf5_filename = "/home/yfz/project/OpenFUSIONToolkit-260316/app/ThinCurr-Python-example-27/tokamak_mesh_EHL2_holes.h5"

try:
    with h5py.File(hdf5_filename, 'r') as f:
        # 检查所需的数据集是否存在
        if 'mesh/R' not in f or 'mesh/LC' not in f or 'mesh/REG' not in f:
            print(f"错误: HDF5 文件 '{hdf5_filename}' 中缺少必要的数据集 ('mesh/R', 'mesh/LC', 'mesh/REG')。")
            exit()
            
        # 读取节点坐标
        nodes = f['mesh/R'][:]
        
        # 读取单元连接 (1-based index)
        elements_1based = f['mesh/LC'][:]
        
        # 读取区域标记
        regions = f['mesh/REG'][:].flatten() # flatten to a 1D array

        num_nodesets = 0
        nodesets = []
        # 检查是否存在 mesh/NUM_NODESETS
        if 'mesh/NUM_NODESETS' not in f:
            print(f"警告: HDF5 文件 '{hdf5_filename}' 中缺少 mesh/NUM_NODESETS 数据集。")
        else:
            # 读取NODESETS数量
            num_nodesets = f['mesh/NUM_NODESETS'][0]

            # 从NODESET0001开始读取
            for i in range(1, num_nodesets + 1):
                nodeset_name = f"NODESET{i:04d}"
                if nodeset_name in f['mesh']:
                    # nodeset = f['mesh'][nodeset_name][:]
                    nodesets.append(f['mesh'][nodeset_name][:])
                    # print(f"  - 读取 {nodeset_name}: 包含 {nodeset.shape[0]} 个节点")
                else:
                    print(f"  - 警告: {nodeset_name} 未找到")

        num_sidesets = 0
        sidesets = []
        # 检查是否存在 mesh/NUM_SIDESETS
        if 'mesh/NUM_SIDESETS' not in f:
            print(f"警告: HDF5 文件 '{hdf5_filename}' 中缺少 mesh/NUM_SIDESETS 数据集。")
        else:
            # 读取SIDESETS数量
            num_sidesets = f['mesh/NUM_SIDESETS'][0]

            # 从SIDSET0001开始读取
            for i in range(1, num_sidesets + 1):
                sideset_name = f"SIDSET{i:04d}"
                if sideset_name in f['mesh']:
                    # sideset = f['mesh'][sideset_name][:]
                    sidesets.append(f['mesh'][sideset_name][:])
                    # print(f"  - 读取 {sideset_name}: 包含 {sideset.shape[0]} 个单元")
                else:
                    print(f"  - 警告: {sideset_name} 未找到")

except FileNotFoundError:
    print(f"错误: 文件 '{hdf5_filename}' 未找到。请先运行 Gmsh 脚本生成该文件。")
    exit()

print(f"成功从 '{hdf5_filename}' 读取数据:")
print(f"  - 节点数量: {nodes.shape[0]}")
print(f"  - 单元数量: {elements_1based.shape[0]}")
print(f"  - 发现的区域标记: {np.unique(regions)}")
print(f"  - 发现的 NODESETS 数量: {num_nodesets}")
print(f"  - 发现的 SIDESETS 数量: {num_sidesets}")

# print(nodesets)


# --- 2. 准备 PyVista 需要的数据格式 ---

# PyVista 需要 0-based 索引，并且单元格式为 [n_points, p0, p1, p2, ...]
# 因为我们所有的单元都是三角形，所以 n_points 总是 3

# 将 1-based 索引转换为 0-based
elements_0based = elements_1based - 1

# 创建 PyVista 格式的单元数组
# 首先，创建一个 (n_elements, 1) 的数组，其中填充了 3
padding = np.full((elements_0based.shape[0], 1), 3, dtype=elements_0based.dtype)

# 将 '3' 列与 0-based 节点索引水平堆叠
# 结果形状为 (n_elements, 4)，例如 [[3, 0, 1, 2], [3, 1, 3, 2], ...]
cells = np.hstack((padding, elements_0based))

# --- 3. 创建 PyVista 网格对象并进行可视化 ---

# 使用节点和单元数据创建 PolyData 对象 (适用于曲面网格)
mesh = pv.PolyData(nodes, faces=cells)

# 将区域标记作为 "cell data" 附加到网格上
mesh.cell_data['Region'] = regions

# 创建一个绘图器
plotter = pv.Plotter()
plotter.set_background('lightgrey') # 设置背景色

# 添加网格到绘图器，并应用可视化样式
plotter.add_mesh(
    mesh,
    scalars='Region',        # 使用 'Region' 数据进行着色
    cmap='viridis',          # 选择一个颜色映射
    show_edges=True,         # 显示网格边
    edge_color='black',      # 边的颜色
    line_width=2,            # 边的线宽
    opacity=1,            # 设置面的半透明度
    scalar_bar_args={'title': 'Region ID'} # 添加颜色图例
)

# 绘制NODESETS

# 设置相机视角
plotter.camera_position = 'xy'
plotter.camera.zoom(1.2)

# 显示交互式窗口
print("\n正在打开 PyVista 可视化窗口...")
plotter.show()
