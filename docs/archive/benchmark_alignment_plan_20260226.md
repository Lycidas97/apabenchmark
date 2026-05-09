# 04_benchmark 路径对齐与修改方案（2026-02-26）

## 1) 目标
以当前 `apabenchmark_final` 既有数据结构为准，保证 `04_benchmark` 的输入/输出路径在逻辑上可与后续分析对接。

## 2) 边界
- 保留 `pf` 数据链路（后续用于 performance）。
- `raw` 的 GT 生成流程暂不并入本项目。
- `raw` 的 BAM 输入来自外部目录：`/path/to/source_raw_bam`。
- `04` 本轮不实际跑全流程，重点是路径映射与结构可连通。

## 3) 路径映射表（当前结构优先）
| 类型 | 上游来源 | 在 `apabenchmark_final` 的统一入口 | 说明 |
|---|---|---|---|
| sim benchmark 输出 | `/path/to/legacy_apabenchmark_root/data/sim_bam_result` | `data/result/benchmark/sim_bam` | 已软链接迁移 |
| raw benchmark 输出 | `/path/to/legacy_apabenchmark_root/data/raw_bam_result` | `data/result/benchmark/raw_bam` | 已软链接迁移 |
| pf benchmark 输出 | 由 `04_benchmark` 运行产出 | `data/result/benchmark/pf_bam` | 保留路径规范，按需后续生成 |
| raw BAM 输入（04） | `/path/to/source_raw_bam` | `data/raw_data/bam/raw_bam` | 新增“重命名映射层” |
| raw GT（analysis 用） | 外部流程导入 | `data/raw_data/raw_bam/{sample}/(bam.expr.tsv,pas.bed)` | 仅数据接入，不纳入生成流程 |
| 旧 raw 输入兼容路径 | - | `data/int_data/bam_to_detect_pas` -> `data/raw_data/bam/raw_bam` | 兼容旧脚本 |

## 4) 已识别冲突与改法
1. `generate_sample_config.sh` 的 raw 输入路径仍指向 `data/int_data/bam_to_detect_pas`（历史路径）。  
   改法：改为统一从 `data/raw_data/bam/raw_bam` 读取，并保留旧路径软链接兼容。

2. raw BAM 命名可能与当前项目样本命名不一致。  
   改法：新增 `snakemake_profile/raw_bam_name_map.tsv`（`canonical_sample -> source_sample`），通过建链脚本做重命名映射。

3. `snakemake_profile/config.yaml` 中 singularity bind 指向旧项目路径。  
   改法：更新为 `/path/to/apabenchmark_final/scripts/04_benchmark`。

4. `pf_bam` 路径在 `rule all` 中被要求，但当前并未提前迁移历史结果。  
   改法：保持 `pf_bam` 作为 `04` 标准产出路径，不做“错误链接”到其他目录，避免语义混淆。

## 5) 新增与改动文件
- `scripts/04_benchmark/link_raw_bam_inputs.sh`  
  从外部 raw_bam 建立标准化输入链接，并校验 raw GT 是否齐备。
- `scripts/04_benchmark/snakemake_profile/raw_bam_name_map.tsv`  
  raw BAM 重命名映射表。
- `scripts/04_benchmark/generate_sample_config.sh`  
  raw 输入目录改为 `data/raw_data/bam/raw_bam`，并增强鲁棒性。
- `scripts/04_benchmark/snakemake_profile/config.yaml`  
  修正 singularity bind 路径。

## 6) 显式假设
1. raw GT 的样本标准名应以 `data/raw_data/raw_bam/` 下目录名为准。  
2. 映射表中的 `canonical_sample` 将作为本项目统一样本名。  
3. 若外部 BAM 命名变化，只需改 `raw_bam_name_map.tsv`，不改 Snakefile。
