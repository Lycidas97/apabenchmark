# Integrated PAS to DaPars2 Converter

## 概述

整合版PAS到DaPars2转换器，使用bedtools建立精确的PAS-基因映射关系，支持独立celltype注释。所有功能集成在一个脚本中，采用模块化设计。

## 使用方法

```bash
source activate scmapa_r
Rscript integrated_pas_to_dapars2_converter.R \
    --pas_counts pas_counts.tsv \
    --pas_coords pas.bed \
    --pas_annotation mm10_sim_pas_gn2000_rep1.bed \
    --celltype_info Chromium_human_embryolimb_...expr.tsv \
    -o output.txt \
    --min_exp 0 \
    --min_cells 0 \
    --debug
```

## 必需参数

- `--pas_counts`: PAS表达矩阵文件 (行: 细胞barcode, 列: PAS ID)
- `--pas_coords`: PAS坐标文件 (BED格式: chr start end peak_id score strand)
- `--pas_annotation`: PAS注释文件 (必须包含exon_id列，格式: ENSMUSG:GeneSymbol:chr:start:end:strand:TE)
- `--celltype_info`: 细胞类型注释文件 (第一列: barcode, 第二列: celltype)

## 可选参数

- `-o, --output`: 输出文件名 (默认: integrated_dapars2_output.txt)
- `--min_exp`: 最小表达阈值 (默认: 0)
- `--min_cells`: 最小细胞数阈值 (默认: 0)
- `--debug`: 启用调试输出

## 输出格式

输出文件包含以下列：
- `Gene`: 基因信息 (格式: ENSMUSG|GeneSymbol|chr|strand)
- `fit_value`: 拟合值 (固定为0)
- `Predicted_Proximal_APA`: 预测近端APA位置
- `Loci`: 基因区域 (格式: chr:start-end)
- `{CellType}_long_exp`: 指定细胞类型的下游PAS表达量
- `{CellType}_short_exp`: 指定细胞类型的上游PAS表达量
- `{CellType}_PDUI`: 指定细胞类型的PDUI值

## 技术特点

1. **自动bedtools调用**: 脚本自动检测和调用bedtools进行坐标映射
2. **模块化设计**: 6个独立模块，便于调试和维护
3. **错误处理**: 完整的输入验证和错误报告
4. **临时文件管理**: 自动清理临时文件
5. **调试支持**: 详细的调试输出选项

## 核心模块

1. **基因区域提取**: 从exon_id列解析基因坐标
2. **PAS-基因映射**: 使用bedtools intersect建立精确映射
3. **数据读取**: 整合表达矩阵和细胞类型信息
4. **PAS信息表**: 创建标准化的PAS注释表
5. **基因处理**: 计算PDUI值和上下游PAS
6. **结果保存**: 输出DaPars2格式文件

## 测试结果

使用测试数据：
- 输入基因数: 1880个
- 处理基因数: 1225个
- 细胞类型: 2个 (T1, T2)
- 输出格式: 完全兼容DaPars2

## 环境要求

- R环境 (scmapa_r conda环境)
- bedtools (自动检测)
- 基础R包 (无需额外安装)

## 注意事项

1. 确保所有输入文件使用相同的坐标系统
2. PAS注释文件必须包含exon_id列
3. 细胞barcode在不同文件间需要保持一致
4. 建议先用小规模数据测试