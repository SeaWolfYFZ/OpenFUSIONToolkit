#!/usr/bin/env python
# coding: utf-8

# In[ ]:


print("Hello World")


# In[ ]:


import numpy as np
import pyvista
import struct


# In[ ]:


from OpenFUSIONToolkit import OFT_env
from OpenFUSIONToolkit.ThinCurr import ThinCurr


# In[ ]:


myOFT = OFT_env(nthreads=16)  # Use 1 thread for visualization
tw_model = ThinCurr(myOFT)


# In[ ]:


restart_dir = "."
hdf5_file_path = restart_dir + "/tokamak_mesh_holes_16.h5"
xml_file_path = restart_dir + "/oft_in.xml"

tw_model.setup_model(mesh_file=hdf5_file_path, xml_filename=xml_file_path)
tw_model.setup_io(basepath=restart_dir)


# In[ ]:


nsteps = 200
tw_model.plot_td(nsteps,compute_B=True,plot_freq=10)
plot_data = tw_model.build_XDMF()


# In[ ]:


pyvista.OFF_SCREEN = True  # 关键：启用离屏渲染
pyvista.start_xvfb()       # 启动虚拟显示（可选，部分场景需要）
pyvista.set_jupyter_backend('static') # Comment to enable interactive PyVista plots

cam_position = [
    (-6, 6, 5),  # 相机位置坐标
    (0.0, 0.0, 0.0),  # 焦点坐标
    (0.0, 0.0, 1.0)  # 定义"向上"方向
]
windowsize = [1920, 1080]


# In[ ]:


grid = plot_data['ThinCurr']['smesh'].get_pyvista_grid()
icoil_grid = plot_data['ThinCurr']['icoils'].get_pyvista_grid()


# In[ ]:


# p = pyvista.Plotter(window_size=windowsize)
# p.camera_position = cam_position
# p.show_axes()
# p.add_mesh(grid, color="white", opacity=0.5, show_edges=False)
# p.add_mesh(grid, style="wireframe", color="black", opacity=1.0, show_edges=True, edge_opacity=1.0)
# p.add_mesh(icoil_grid, color="blue", opacity=1.0, line_width=4)
# # p.add_ruler(pointa=(0, -1.5, 0), pointb=(1.5, -1.5, 0), title='')
# # p.add_box_axes()

# p.show()


# In[ ]:


# pyvista.OFF_SCREEN = True  # 关键：启用离屏渲染
# pyvista.start_xvfb()       # 启动虚拟显示（可选，部分场景需要）
# pyvista.set_jupyter_backend('static') # Comment to enable interactive PyVista plots

# cam_position = [
#     (-6, 6, 5),  # 相机位置坐标
#     (0.0, 0.0, 0.0),  # 焦点坐标
#     (0.0, 0.0, 1.0)  # 定义"向上"方向
# ]
# windowsize = [1920, 540]  # 宽度增加3倍以容纳三个子图

# plot_time = 1.0E-3
# plot_step = 10
# Jfull = plot_data['ThinCurr']['smesh'].get_field('J_v',timestep = plot_step)
# Bfull = plot_data['ThinCurr']['smesh'].get_field('B_v',timestep = plot_step)
# Ffull = np.zeros((np.shape(Jfull)))
# print(np.shape(Ffull))
# for i in range(np.shape(Jfull)[0]):
#     Ffull[i,0] = Jfull[i,1] * Bfull[i,2] - Jfull[i,2] * Bfull[i,1]
#     Ffull[i,1] = Jfull[i,2] * Bfull[i,0] - Jfull[i,0] * Bfull[i,2]
#     Ffull[i,2] = Jfull[i,0] * Bfull[i,1] - Jfull[i,1] * Bfull[i,0] 

# J_scale = 0.25 / 4e5
# B_scale = 0.5 / 3
# F_scale = 0.25 / 2.5e5

# # 创建3个子图的Plotter
# p = pyvista.Plotter(shape=(1, 3), window_size=windowsize)

# # 第一个子图：J
# grid["vectors"] = Jfull
# grid.set_active_vectors("vectors")
# p.subplot(0, 0)
# arrows = grid.glyph(scale="vectors", orient="vectors", factor=J_scale)
# p.camera_position = cam_position
# # p.show_axes()
# p.add_mesh(grid, color="white", opacity=0.25, show_edges=False)
# p.add_mesh(arrows, cmap="turbo", scalar_bar_args={'title': "J", "vertical": True, "position_x": colorbar_position_x})

# # 第二个子图：B
# grid["vectors"] = Bfull
# grid.set_active_vectors("vectors")
# p.subplot(0, 1)
# arrows = grid.glyph(scale="vectors", orient="vectors", factor=B_scale)
# p.camera_position = cam_position
# # p.show_axes()
# p.add_mesh(grid, color="white", opacity=0.25, show_edges=False)
# p.add_mesh(arrows, cmap="turbo", scalar_bar_args={'title': "B", "vertical": True, "position_x": colorbar_position_x})

# # 第三个子图：F=J×B
# grid["vectors"] = Ffull
# grid.set_active_vectors("vectors")
# p.subplot(0, 2)
# arrows = grid.glyph(scale="vectors", orient="vectors", factor=F_scale)
# p.camera_position = cam_position
# # p.show_axes()
# p.add_mesh(grid, color="white", opacity=0.25, show_edges=False)
# p.add_mesh(arrows, cmap="turbo", scalar_bar_args={'title': "F=J×B", "vertical": True, "position_x": colorbar_position_x})

# # 显示所有子图
# p.show()


# In[ ]:


import subprocess
import os


# In[ ]:


import pyvista as pv
import numpy as np
import os

# ==========================================
# 1. 读取 C++ 导出的环面磁场数据 (torus_B.bin)
# ==========================================
def read_torus_b(path):
    with open(path, "rb") as f:
        magic = f.read(8)
        if magic != b"OFTTORUS":
            raise ValueError(f"Bad magic in {path}: {magic!r}")

        version = struct.unpack("<i", f.read(4))[0]
        if version != 1:
            raise ValueError(f"Unsupported torus_B.bin version: {version}")

        nphi = struct.unpack("<i", f.read(4))[0]
        ntheta = struct.unpack("<i", f.read(4))[0]
        npts = struct.unpack("<i", f.read(4))[0]
        nsaves = struct.unpack("<i", f.read(4))[0]
        save_stride_steps = struct.unpack("<i", f.read(4))[0]
        R0 = struct.unpack("<d", f.read(8))[0]
        a = struct.unpack("<d", f.read(8))[0]
        dt = struct.unpack("<d", f.read(8))[0]

        pts = np.fromfile(f, dtype="<f8", count=3 * npts).reshape((npts, 3))

        times = np.empty((nsaves,), dtype=np.float64)
        B = np.empty((nsaves, npts, 3), dtype=np.float64)
        for i in range(nsaves):
            times[i] = struct.unpack("<d", f.read(8))[0]
            B[i, :, :] = np.fromfile(f, dtype="<f8", count=3 * npts).reshape((npts, 3))

    return {
        "nphi": nphi,
        "ntheta": ntheta,
        "npts": npts,
        "nsaves": nsaves,
        "save_stride_steps": save_stride_steps,
        "R0": R0,
        "a": a,
        "dt": dt,
        "pts": pts,
        "times": times,
        "B": B,
    }


# ==========================================
# 2. 参数设置与预计算
# ==========================================
# 设置离屏渲染
pv.OFF_SCREEN = True
# pv.start_xvfb() # 如果在无头服务器上运行请取消注释
# pv.set_jupyter_backend('static') # 如果在Jupyter中运行

cam_position = [(-6, 6, 5), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
windowsize = [1920, 1080] # 调整为更适合 2x2 的比例

# 时间参数
torus_bin_file = os.path.join(restart_dir, "torus_B.bin")
torus_data = read_torus_b(torus_bin_file)
torus_pts = torus_data["pts"]
torus_times = torus_data["times"]
torus_B = torus_data["B"]
torus_R0 = torus_data["R0"]

dt = float(torus_data["dt"])
plot_freq = int(torus_data["save_stride_steps"])  # should match C++ export stride and OFT plot_td plot_freq
plot_dt = dt * plot_freq
start_step = 0
end_step = 40
end_step = min(end_step, int(torus_data["nsaves"]) - 1)

# --- 手动设置统一的标量范围 (User Requirement 2) ---
# 请根据你的实际物理量级修改这些值
J_clim = [0, 6.0e5]   # 例如：J 的范围 0 到 2 MA/m^2
B_clim = [0, 3.5]     # 例如：B 的范围 0 到 3 T
F_clim = [0, 8.0e5]   # 例如：F 的范围
# ------------------------------------------------

# 缩放因子 (用于箭头大小)
J_scale = 0.25 / 4e5
B_scale = 0.5 / 3
F_scale = 0.25 / 2.5e5

colorbar_position_x = 0.01

fig_dir = "./figs"
os.makedirs(fig_dir, exist_ok=True)

# Build a torus surface mesh connectivity once (2 triangles per param quad).
nphi = int(torus_data["nphi"])
ntheta = int(torus_data["ntheta"])
faces = []
for iphi in range(nphi):
    iphi2 = (iphi + 1) % nphi
    for itheta in range(ntheta):
        itheta2 = (itheta + 1) % ntheta
        p00 = iphi * ntheta + itheta
        p10 = iphi2 * ntheta + itheta
        p01 = iphi * ntheta + itheta2
        p11 = iphi2 * ntheta + itheta2
        faces.extend([3, p00, p10, p11])
        faces.extend([3, p00, p11, p01])
faces = np.array(faces, dtype=np.int64)
torus_mesh = pv.PolyData(torus_pts, faces)

# Precompute outward unit normals for the torus surface:
# n = (p - centerline(phi)) / ||p - centerline(phi)||
phi = np.arctan2(torus_pts[:, 1], torus_pts[:, 0])
centerline = np.stack([torus_R0 * np.cos(phi), torus_R0 * np.sin(phi), np.zeros_like(phi)], axis=1)
nvec = torus_pts - centerline
nvec /= np.linalg.norm(nvec, axis=1)[:, None]

# Global color limits for Bn across all saved times, for consistent scaling.
Bn_all = np.einsum("tnc,nc->tn", torus_B, nvec)
Bn_clim = float(np.max(np.abs(Bn_all)))

# ==========================================
# 3. 循环绘制
# ==========================================
for plot_step in range(start_step, end_step + 1):
    print(f"正在处理时间步 {plot_step}/{end_step}...")

    # 获取数据
    try:
        Jfull = plot_data['ThinCurr']['smesh'].get_field('J_v', timestep=plot_step)
        Bfull = plot_data['ThinCurr']['smesh'].get_field('B_v', timestep=plot_step)
    except Exception as e:
        print(f"获取时间步 {plot_step} 的数据失败: {e}")
        continue

    # 1. 计算 F 向量 (F = J × B)
    Ffull = np.zeros_like(Jfull)
    Ffull[:, 0] = Jfull[:, 1] * Bfull[:, 2] - Jfull[:, 2] * Bfull[:, 1]
    Ffull[:, 1] = Jfull[:, 2] * Bfull[:, 0] - Jfull[:, 0] * Bfull[:, 2]
    Ffull[:, 2] = Jfull[:, 0] * Bfull[:, 1] - Jfull[:, 1] * Bfull[:, 0]

    # 2. [新增] 计算力的大小（模长）用于热力图
    F_mag = np.linalg.norm(Ffull, axis=1)
    J_norm = np.linalg.norm(Jfull, axis=1)

    current_time = float(torus_times[plot_step])

    # 创建 2x2 布局
    p = pv.Plotter(shape=(2, 2), window_size=windowsize, off_screen=True)

    # --- 子图 1 (0,0): J (保持矢量箭头) ---
    p.subplot(0, 0)
    grid["vectors"] = Jfull
    grid.set_active_vectors("vectors")
    arrows_J = grid.glyph(scale="vectors", orient="vectors", factor=J_scale)
    max_J = np.max(J_norm)
    p.add_mesh(grid, color="white", opacity=0.25, show_edges=False)
    p.add_mesh(arrows_J, cmap="turbo", clim=J_clim,
               scalar_bar_args={'title': "J(A/m^2)", "vertical": True, "n_labels": 5, "position_x": colorbar_position_x})
    p.add_text(f"Current Density J (Max: {max_J:.2e})\nt={current_time*1000:.2f} ms", font_size=10)
    p.camera_position = cam_position

    # --- 子图 2 (0,1): B (保持矢量箭头) ---
    p.subplot(0, 1)
    grid["vectors"] = Bfull
    grid.set_active_vectors("vectors")
    arrows_B = grid.glyph(scale="vectors", orient="vectors", factor=B_scale)
    p.add_mesh(grid, color="white", opacity=0.25, show_edges=False)
    p.add_mesh(arrows_B, cmap="turbo", clim=B_clim,
               scalar_bar_args={'title': "B(T)", "vertical": True, "n_labels": 5, "position_x": colorbar_position_x})
    p.add_text("Magnetic Field B", font_size=10)
    p.camera_position = cam_position

    # --- 子图 3 (1,0): F (修改为：力的大小热力图) ---
    p.subplot(1, 0)

    # 将计算好的模长赋值给网格
    grid["F_Magnitude"] = F_mag

    max_f = np.max(F_mag)
    print(f"  -> Time: {current_time*1000:.2f}ms, Max Force: {max_f:.4e}")

    # 直接绘制网格表面，使用 scalars 参数进行着色
    # cmap="magma" 或 "hot" 非常适合表示“力/能量”的大小
    # clim=F_clim 需要确保您在外部定义了 F_clim (例如 [0, max_force])，否则去掉该参数让其自动缩放
    p.add_mesh(grid, 
               scalars="F_Magnitude", 
               cmap="turbo",       # 改用 turbo，0值为深蓝，最大值为红
               clim=F_clim,        # 确保 F_clim 的范围与 print 出来的 Max Force 数量级匹配
               lighting=False,     # <--- 关键修改：关闭光照
               show_edges=False,   
               scalar_bar_args={'title': "F(N/m^2)", 
                                "vertical": True, 
                                "n_labels": 5, 
                                "position_x": colorbar_position_x,
                                "color": "black"} # 如果背景是白色的，字设为黑；背景黑则设为白
               )

    p.add_text(f"Lorentz force |F| (Max: {max_f:.2e})", font_size=10)
    p.camera_position = cam_position

    # --- 子图 4 (1,1): 环面上的 B_n 可视化 (替换掉不正确的输入波形图) ---
    p.subplot(1, 1)
    Bn = Bn_all[plot_step, :]
    torus_mesh.point_data["B_n"] = Bn
    p.add_mesh(
        torus_mesh,
        scalars="B_n",
        cmap="coolwarm",
        clim=[-Bn_clim, Bn_clim],
        lighting=False,
        show_edges=False,
        scalar_bar_args={
            "title": "Torus B_n (T)",
            "vertical": True,
            "n_labels": 5,
            "position_x": colorbar_position_x,
            "color": "black",
        },
    )
    p.add_text(f"Torus B_n on surface (t={current_time*1000:.2f} ms)", font_size=10)
    p.camera_position = cam_position

    # 保存图片
    filename = os.path.join(fig_dir, f"frame_{plot_step:04d}.png")
    p.screenshot(filename)
    p.close()

    print(f"已保存: {filename}")


print(f"所有图片已保存到 {fig_dir} 目录")


# In[ ]:


# 使用ffmpeg将图片编码为视频
print("\n正在使用ffmpeg生成视频...")
# 创建保存视频的目录
video_dir = "./anime"
os.makedirs(video_dir, exist_ok=True)
# 视频文件名
video_path = os.path.join(video_dir, "simulation.mp4")
# ffmpeg命令
# -r 30: 设置帧率为30fps
# -i: 输入文件模式
# -c:v libx264: 使用H.264编码
# -pix_fmt yuv420p: 确保兼容性
# -vf "fps=30": 确保输出帧率
ffmpeg_cmd = [
    'ffmpeg',
    '-y',  # 覆盖输出文件而不询问
    '-framerate', '10',  # 输入帧率
    '-i', os.path.join(fig_dir, 'frame_%04d.png'),  # 输入文件模式
    '-c:v', 'libx264',  # 视频编码器
    '-pix_fmt', 'yuv420p',  # 像素格式
    '-vf', 'fps=30',  # 输出帧率
    '-crf', '18',  # 质量参数（0-51，值越小质量越高）
    '-preset', 'slow',  # 编码速度预设
    video_path
]
try:
    # 执行ffmpeg命令
    print("执行命令:", ' '.join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True, check=True)
    print("视频生成成功!")
    print(f"视频已保存到: {video_path}")

    # 如果需要，可以添加更多视频编码选项
    # 例如，生成GIF动画：
    # gif_path = os.path.join(video_dir, "simulation.gif")
    # gif_cmd = [
    #     'ffmpeg',
    #     '-y',
    #     '-framerate', '30',
    #     '-i', os.path.join(fig_dir, 'frame_%04d.png'),
    #     '-vf', 'fps=15,scale=960:-1:flags=lanczos',
    #     '-c:v', 'gif',
    #     gif_path
    # ]
    # subprocess.run(gif_cmd, check=True)
    # print(f"GIF动画已保存到: {gif_path}")

except subprocess.CalledProcessError as e:
    print(f"ffmpeg执行失败: {e}")
    print(f"标准错误输出:\n{e.stderr}")
except FileNotFoundError:
    print("错误: 未找到ffmpeg。请确保已安装ffmpeg并将其添加到系统PATH中。")
    print("在Ubuntu上可以使用: sudo apt install ffmpeg")
    print("在macOS上可以使用: brew install ffmpeg")
    print("在Windows上可以从官网下载: https://ffmpeg.org/download.html")
