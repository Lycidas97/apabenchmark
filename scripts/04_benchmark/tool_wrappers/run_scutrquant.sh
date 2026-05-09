#!/bin/bash

# 帮助文档
function display_help() {
    echo "Usage: $0 [options]"
    echo
    echo "   -s, --sample_file FILE     样本文件路径"
    echo "   -t, --tmp_dir DIR          临时目录路径"
    echo "   -c, --targets_config FILE  目标配置文件"
    echo "   -r, --target STRING        分析目标"
    echo "   -h, --help                 显示帮助信息"
    echo
    exit 1
}

# 参数解析
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--sample_file) sample_file="$2"; shift ;;
        -t|--tmp_dir) tmp_dir="$2"; shift ;;
        -c|--targets_config) targets_config="$2"; shift ;;
        -r|--target) target="$2"; shift ;;
        -h|--help) display_help ;;
        *) echo "Unknown parameter: $1"; display_help ;;
    esac
    shift
done

# 参数检测
if [ -z "$sample_file" ] || [ -z "$tmp_dir" ] || [ -z "$targets_config" ] || [ -z "$target" ]; then
    echo "Error: Missing required parameter (sample_file, tmp_dir, targets_config, target)"
    display_help
fi

echo 