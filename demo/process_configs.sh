#!/bin/bash

# ==========================================================
# Shell Script: process_configs.sh
# Usage: ./process_configs.sh <input_dir> <output_dir>
# Description: Reads *.cfg files from <input_dir>, parses parameters, 
#              and generates parameter files in <output_dir>/<cfg_name>/par/.
#              It calls 'generate_ViewAngleList.py' to create the angle list.
# ==========================================================

# 检查参数数量
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "  <input_dir>: 包含所有配置文件的文件夹 (例如: ./configs)"
    echo "  <output_dir>: 存放所有输出结果的根文件夹 (例如: ./output)"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
genViewList="./generate_ViewAngleList.py"

# 检查输入文件夹是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' not found."
    exit 1
fi

# 检查 Python 脚本是否存在（关键步骤，确保能调用）
if [ ! -f "generate_ViewAngleList.py" ]; then
    echo "Error: Python script 'generate_ViewAngleList.py' not found in the current directory."
    echo "Please ensure your Python script is created and placed here."
    exit 1
fi

# 创建输出根文件夹 (如果不存在)
mkdir -p "$OUTPUT_DIR"
echo "Output root directory created: $OUTPUT_DIR"
echo "--------------------------------------------------"

# --- 函数：读取配置文件中的值 ---
# 参数1: 配置文件路径
# 参数2: 变量名 (e.g., NViews)
get_config_value() {
    # 使用 grep 查找变量行，并用 awk 提取等号后的值，再去除首尾空白
    # 查找以变量名开头，后面跟等号的行，忽略注释
    grep "^[[:space:]]*$2[[:space:]]*=" "$1" | awk -F'=' '{print $2}' | tr -d '[:space:]'
}

# --- 循环处理输入文件夹中所有的 .cfg 文件 ---
for cfg_file in "$INPUT_DIR"/*.cfg; do
    
    # 防止在没有匹配文件时，循环处理 "*.cfg" 字符串本身
    if [ ! -f "$cfg_file" ]; then
        if [[ "$cfg_file" == "$INPUT_DIR/*.cfg" ]]; then
            echo "Warning: No .cfg files found in '$INPUT_DIR'. Exiting loop."
        fi
        continue
    fi
    
    # 获取文件名，不包含路径和扩展名 (e.g., r2a)
    base_name=$(basename "$cfg_file" .cfg)
    
    # 定义当前的输出文件夹 (e.g., ./output/r2a)
    current_output_dir="${OUTPUT_DIR}/${base_name}"
    par_dir="${current_output_dir}/par"
    
    echo ">>> Processing $cfg_file -> $current_output_dir <<<"

    # 1. 解析配置文件
    NViews=$(get_config_value "$cfg_file" "NViews")
    NChannels=$(get_config_value "$cfg_file" "NChannels")
    DeltaChannel=$(get_config_value "$cfg_file" "DeltaChannel")
    CenterOffset=$(get_config_value "$cfg_file" "CenterOffset")
    Nx=$(get_config_value "$cfg_file" "Nx")
    Ny=$(get_config_value "$cfg_file" "Ny")
    DeltaPix=$(get_config_value "$cfg_file" "DeltaPix")

    # 检查关键参数是否被成功读取
    if [ -z "$NViews" ] || [ -z "$Nx" ] || [ -z "$DeltaPix" ] || [ -z "$NChannels" ] || [ -z "$DeltaChannel" ] || [ -z "$CenterOffset" ]; then
        echo "Error: Could not read all required parameters from $cfg_file. Check your config file format. Skipping."
        continue
    fi

    # 2. 创建目录结构
    mkdir -p "${par_dir}"
    echo "  - Created directory structure: ${par_dir}"
    
    # 3. 生成 r2a.imgparams
    imgparams_file="${par_dir}/${base_name}.imgparams"
    
    cat << EOF > "$imgparams_file"
Nx: ${Nx} 
Ny: ${Ny}
Nz: 1
FirstSliceNumber: 1
Deltaxy: ${DeltaPix}
DeltaZ: 0.1
ROIRadius: 250
EOF
    
    echo "  - Generated ${base_name}.imgparams"

    # 4. 生成 r2a.sinoparams
    sinoparams_file="${par_dir}/${base_name}.sinoparams"
    
    cat << EOF > "$sinoparams_file"
Geometry : parallel
NChannels: ${NChannels} 
NViews: ${NViews} 
NSlices: 1
DeltaChannel: ${DeltaChannel}
CenterOffset: ${CenterOffset}
DeltaSlice: 0.9765625
FirstSliceNumber: 1
ViewAngleList: ./ViewAngleList.txt
EOF

    echo "  - Generated ${base_name}.sinoparams"

    # 5. 调用用户提供的 Python 脚本生成 ViewAngleList.txt
    viewangle_file="${par_dir}/ViewAngleList.txt"
    
    # **核心调用：将 NViews 作为参数传给 Python 脚本，并重定向输出**
    python $genViewList "${NViews}" > "$viewangle_file"
    
    # 检查 Python 脚本是否成功执行
    if [ $? -eq 0 ]; then
        echo "  - Generated ViewAngleList.txt using NViews=${NViews}"
    else
        echo "  - ERROR: Python script 'generate_ViewAngleList.py' failed for NViews=${NViews}."
        echo "  - Please check your Python script for errors."
    fi

    echo "--- Finished $base_name ---"
    echo ""

done

echo "=================================================="
echo "All configuration files processed and output to: $OUTPUT_DIR"
echo "=================================================="
