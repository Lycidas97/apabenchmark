# Raw Dataset Inventory

This manifest documents the public data resources used to prepare the Stage 02
peak-extraction BAM inputs. The BAM payloads themselves are not committed to
git. Stage 02 expects locally prepared files at:

- `data/raw_data/bam_for_peak_extraction/{sample}.bam`
- `data/raw_data/bam_for_peak_extraction/{sample}.bam.bai`

The inventory is dataset/resource-level rather than per-BAM. The project Table
S1 workbook records 32 resource/accession entries: 23 empirical
simulation-design studies, 2 paired long-read validation studies, and additional
public 10x/STOmics demo resources. The 289 Stage 02 samples are summarized by
sample prefix below from `scripts/02_peak_analysis/snakemake_config.yaml`.

## Public Resource Inventory

| Study ID | Resource/accession | Repository | Source URL | Short label | Protocol/assay | Species | Tissue/source | Benchmark role | Release status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SIM-01 | GSE121891 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121891` | chromium_mouse_ob | 10x Chromium | Mouse | Olfactory bulb | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-02 | GSE184708 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184708` | chromium_mouse_testicle | 10x Chromium | Mouse | Testicle/gonad | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-03 | GSE185671 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185671` | chromium_mouse_retina | 10x Chromium | Mouse | Retina | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-04 | GSE120716 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120716` | chromium_human_prostate | 10x Chromium | Human | Prostate | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-05 | GSE222322 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222322` | visium_human_spine | 10x Visium | Human | Spinal cord/spine | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-06 | GSE117011 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117011` | dropseq_mouse_heart | Drop-seq | Mouse | Heart valve/heart | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-07 | GSE112393 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112393` | dropseq_mouse_testicle | Drop-seq | Mouse | Testicle | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-08 | GSE148829 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148829` | dropseq_mouse_epicell | Drop-seq | Mouse | Epithelial cells | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-09 | GSE124872 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124872` | dropseq_mouse_lung | Drop-seq | Mouse | Lung | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-10 | GSE97942 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942` | dropseq_human_brain | snDrop-seq / Drop-seq | Human | Brain | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-11 | GSE129611 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129611` | human_hairfollicle | Drop-seq / 10x Chromium | Human | Hair follicle | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-12 | GSE153562 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153562` | microwell_atlas | Microwell-seq | Mouse | Brain/testicle atlas subsets | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-13 | GSE134355 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134355` | microwell_human_embryo | Microwell-seq | Human | Embryo | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-14 | GSE120374 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120374` | spatialtranscriptome_mouse_spine | Spatial Transcriptomics | Mouse | Spine | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-15 | GSE169021 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169021` | slideseq_mouse_ob | Slide-seq V2 | Mouse | Olfactory bulb | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-16 | PRJNA668433 | NCBI BioProject | `https://www.ncbi.nlm.nih.gov/bioproject/PRJNA668433` | slideseq_mouse_testicle | Slide-seq V2 | Mouse | Testicle | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-17 | PRJNA821403 | NCBI BioProject / CNGBdb | `https://www.ncbi.nlm.nih.gov/bioproject/PRJNA821403` | visium_human_liver | 10x Visium | Human | Liver | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-18 | E-MTAB-8813 | ArrayExpress / BioStudies | `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8813` | chromium_human_embryo | 10x Chromium | Human | Embryo limb | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-19 | SCP815 | Broad Single Cell Portal | `https://singlecell.broadinstitute.org/single_cell/study/SCP815` | slideseq_mouse_dataset | Slide-seq V2 | Mouse | Hippocampus / embryo / brain subsets | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-20 | STDS0000058 | STOmicsDB | `https://db.cngb.org/stomics/project/STDS0000058` | stereoseq_mosta | Stereo-seq | Mouse | MOSTA embryo atlas | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-21 | Stereo-seq mouse brain demo C03937C4 | STOmics demo | `https://en.stomics.tech/col1201/index` | stereoseq_mouse_brain_C03937C4 | Stereo-seq | Mouse | Brain | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-22 | Stereo-seq mouse brain demo C04042E3 | STOmics demo | `https://en.stomics.tech/col1201/index` | stereoseq_mouse_brain_C04042E3 | Stereo-seq | Mouse | Brain | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-23 | Stereo-seq mouse kidney demo | STOmics demo | `https://en.stomics.tech/col1201/index` | stereoseq_mouse_kidney | Stereo-seq | Mouse | Kidney | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-24 | Stereo-seq mouse ovary demo | STOmics demo | `https://en.stomics.tech/col1201/index` | stereoseq_mouse_ovary | Stereo-seq | Mouse | Ovary | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-25 | Stereo-seq mouse tongue demo | STOmics demo | `https://en.stomics.tech/col1201/index` | stereoseq_mouse_tongue | Stereo-seq | Mouse | Tongue | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-26 | 10x Visium mouse brain sagittal anterior | 10x Genomics | `https://www.10xgenomics.com/datasets` | visium_mouse_brain_sagittal_anterior | 10x Visium | Mouse | Brain | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-27 | 10x Visium mouse brain sagittal posterior | 10x Genomics | `https://www.10xgenomics.com/datasets` | visium_mouse_brain_sagittal_posterior | 10x Visium | Mouse | Brain | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-28 | 10x Visium mouse olfactory bulb | 10x Genomics | `https://www.10xgenomics.com/datasets` | visium_mouse_olfactory_bulb | 10x Visium | Mouse | Olfactory bulb | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-29 | 10x Visium mouse kidney | 10x Genomics | `https://www.10xgenomics.com/datasets` | visium_mouse_kidney | 10x Visium | Mouse | Kidney | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| SIM-30 | Spatial Transcriptomics mouse olfactory bulb | Spatial Research | `https://www.spatialresearch.org/resources-published-datasets/` | st_mouse_olfactory_bulb | Spatial Transcriptomics | Mouse | Olfactory bulb | Empirical peak characterization and simulation design | Public source documented; prepared BAM payload not included |
| VAL-01 | PRJNA599962 | NCBI BioProject | `https://www.ncbi.nlm.nih.gov/bioproject/PRJNA599962` | pbmc_10x_r2c2 | 10x Chromium with R2C2 long-read validation | Human | Peripheral blood mononuclear cells | Paired short-read/long-read validation | Public source documented; rawbenchmark payload deferred |
| VAL-02 | GSE158450 | GEO | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158450` | mouse_brain_lr_10x_singlecell | 10x Chromium with PacBio long-read validation | Mouse | P7 hippocampus / prefrontal cortex | Paired short-read/long-read validation | Public source documented; rawbenchmark payload deferred |

## Stage 02 Sample Prefix Summary

The Stage 02 peak-analysis configuration currently contains 289 empirical
samples. Counts below are grouped by the first three underscore-separated
tokens of the sample names in `scripts/02_peak_analysis/snakemake_config.yaml`.

| Sample prefix | Count |
| --- | ---: |
| Dropseq_human_brain | 86 |
| Dropseq_mouse_testicle | 41 |
| Visium_human_spinalcord | 18 |
| Chromium_human_prostate | 16 |
| Chromium_mouse_retina | 13 |
| Chromium_human_embryolimb | 12 |
| Microwell_human_embryo | 12 |
| Slideseq_mouse_olfactorybulb | 12 |
| SpatialTranscriptomics_mouse_olfactorybulb | 11 |
| Stereoseq_mouse_embryo | 10 |
| Dropseq_mouse_lung | 8 |
| Dropseq_mouse_epithelial | 7 |
| Slideseq_mouse_testicle | 6 |
| Slideseq_mouse_hippocampus | 5 |
| Chromium_mouse_olfactorybulb | 4 |
| Dropseq_human_hairfollicle | 3 |
| Chromium_human_hairfollicle | 2 |
| Chromium_mouse_testicle | 2 |
| Dropseq_mouse_heart | 2 |
| Microwell_mouse_brain | 2 |
| Microwell_mouse_testicle | 2 |
| SpatialTranscriptomics_mouse_spine | 2 |
| Stereoseq_mouse_olfactorybulb | 2 |
| Visium_mouse_brain | 2 |
| Slideseq_mouse_E15-brain | 1 |
| Slideseq_mouse_embryo | 1 |
| Stereoseq_mouse_brain | 1 |
| Stereoseq_mouse_kidney | 1 |
| Stereoseq_mouse_ovary | 1 |
| Stereoseq_mouse_tongue | 1 |
| Visium_human_liver | 1 |
| Visium_mouse_kidney | 1 |
| Visium_mouse_olfactorybulb | 1 |

## Release Boundaries

- Prepared BAM and BAI files are external runtime inputs and are not included in
  the initial git release.
- Per-BAM checksums are not recorded for the initial release because the
  prepared BAM payload is not being published.
- Raw benchmark inputs and paired long-read benchmark payloads remain deferred;
  they are documented here only as source studies.
- If a future release publishes prepared BAM payloads, add an external download
  location and per-file checksums for the exact released BAM/BAI files.
