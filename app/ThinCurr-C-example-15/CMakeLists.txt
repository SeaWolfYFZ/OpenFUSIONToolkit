cmake_minimum_required(VERSION 3.10)
project(ThinCurrExample)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 输出目录设置
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# OFT库配置
set(OFT_INSTALL_DIR "/home/yfz/project/OpenFUSIONToolkit-251222/install_release")
set(OFT_LIB "${OFT_INSTALL_DIR}/bin/liboftpy.so")

# 检查OFT库是否存在
if(NOT EXISTS ${OFT_LIB})
    message(WARNING "OFT library not found at ${OFT_LIB}")
endif()

# 包含目录
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# 可执行文件
add_executable(tokamak-test
    src/main.cpp
)

# 链接库
target_link_libraries(tokamak-test ${OFT_LIB})

# 设置运行时库路径
set_target_properties(tokamak-test PROPERTIES
    INSTALL_RPATH "${OFT_INSTALL_DIR}/bin"
    BUILD_WITH_INSTALL_RPATH TRUE
    SKIP_BUILD_RPATH FALSE
)

# 安装目标（可选）
install(TARGETS tokamak-test
    RUNTIME DESTINATION bin
)