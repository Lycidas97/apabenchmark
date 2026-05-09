# GitHub Release Checklist

This checklist tracks the remaining work before a formal public GitHub release.

## Metadata

| Status | Item |
| --- | --- |
| Ready | MIT `LICENSE` is present. |
| Ready | Repository-level pre-release `CITATION.cff` points to `https://github.com/Lycidas97/apabenchmark`. |
| Deferred | Replace repository-level citation with manuscript citation when DOI or preprint details are available. |

## Data and SIF Inputs

| Status | Item |
| --- | --- |
| Ready | Public annotation and mouse FASTA sources are listed in `docs/external_inputs_manifest.md`; PolyASite and PolyA_DB download/extracted-file checksums are recorded. |
| Ready | Stage 02 peak-extraction BAM sources are documented at dataset level in `docs/raw_dataset_inventory.md`; prepared BAM/BAI payloads are external and not included in the initial release. |
| Ready | Tool-reference annotations are documented per file and can be rebuilt from public GENCODE inputs with `scripts/04_benchmark/prepare_tool_references.py`; SCAPE BEDs require a SCAPE container or checkout. |
| Ready | SIF download URLs are documented and SHA256 checksums are recorded from the source image files. |
| Deferred | Raw benchmark inputs are not included in the initial public release. |
| Ready | Large runtime data are ignored by git and documented as external inputs. |

## Workflow

| Status | Item |
| --- | --- |
| Ready | Top-level `scripts/run_full_repro.sh` coordinates stages 01-06. |
| Ready | Stage 01 has a script entrypoint; notebooks are references only. |
| Deferred | Run `./scripts/run_full_repro.sh --dry-run` from a shell with conda, Snakemake, external inputs, and Singularity/Apptainer available. |
| Deferred | Run at least one full or staged reproduction pass before tagging a paper-grade release. |

## Documentation

| Status | Item |
| --- | --- |
| Ready | README describes the repository/data split, quick start, main workflow, outputs, validation limits, and citation/license state. |
| Ready | `docs/full_repro_data_flow.md` defines the connected data flow and stage commands. |
| Ready | Stage READMEs point users back to the top-level runner and current data contract. |
| Optional | Add `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, or `SECURITY.md` if the project will accept external issues or pull requests. |

## Repository Hygiene

| Status | Item |
| --- | --- |
| Ready | No generated BAM/SIF/PDF/log/pid/parquet/cache files are tracked. |
| Ready | Active docs and scripts do not contain private workstation paths. |
| Ready | Archive files are labeled as historical and may contain obsolete examples. |

## Validation Commands

Run these before tagging:

```bash
git status --short
rg -n --fixed-strings --file docs/release_blocklist.txt \
  README.md docs scripts CITATION.cff -g '!*.ipynb' -g '!docs/archive/**' \
  -g '!docs/release_blocklist.txt'
git ls-files '*.sh' scripts | rg '\.sh$' | xargs -r bash -n
git ls-files '*.py' scripts docs/archive/legacy_simulation_scripts | rg '\.py$' | xargs -r python3 -m py_compile
ruby -e 'require "yaml"; YAML.load_file("CITATION.cff")'
```
