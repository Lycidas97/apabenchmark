Benchmarking
============

This section describes the steps to run the benchmarking workflow using Snakemake.

For release reproduction, prefer the top-level `scripts/run_full_repro.sh`
entrypoint. Use the commands below for Stage 04 debugging or reruns after
Stages 01-03 have produced simulation inputs.

Prerequisites
-------------

- Make sure you have Snakemake installed. If not, you can install it using:

  conda install -c bioconda -c conda-forge snakemake

- Make sure you have the necessary dependencies installed. Snakemake will automatically create and manage the required environments based on the defined rules.

Steps
-----

1. Download the required SIF image by running the following command:

   bash download_sif.sh

   This script downloads the SIF images needed for the benchmarking workflow.
   Record release checksums in `docs/external_inputs_manifest.md` after
   downloading final public images.

2. Prepare the tool reference annotations by running:

   python3 prepare_tool_references.py \
     --scape-image ./sif/scape.sif

   This downloads public GENCODE GTF/FASTA inputs, writes the GTF and chromosome
   size references, derives the scMAPA-compatible terminal-exon BED files, and
   calls SCAPE's official `prepare` command to build the SCAPE UTR/intron BED
   files. Use `--scape-repo /path/to/SCAPE` instead of `--scape-image` when
   running from a local SCAPE checkout. Use `--skip-scape` only for a partial
   preparation when SCAPE is not available. Recompute
   `docs/external_inputs_manifest.md` checksums whenever these external payloads
   are regenerated or replaced.

3. Prepare raw BAM links (optional but recommended) by running:

   bash link_raw_bam_inputs.sh /path/to/source_raw_bam

   This script maps external raw BAM files into `data/raw_data/bam/raw_bam` using
   `snakemake_profile/raw_bam_name_map.tsv` and creates a legacy-compatible link
   at `data/int_data/bam_to_detect_pas`.

4. Generate a full config (raw+sim+pf) for `Snakefile` by running:

   bash generate_sample_config.sh /path/to/apabenchmark_final \
     ./snakemake_profile/sample.yaml

   This script writes a complete config to `snakemake_profile/sample.yaml`, which is
   consumed by `configfile: "snakemake_profile/sample.yaml"` in `Snakefile`.

   Note: stage 04 now supports mixed `mm10`/`hg38` sim samples by selecting
   genome-matched references from sample names. You can restrict/override via
   `supported_sim_genomes` in `snakemake_profile/sample.yaml` (e.g. `["mm10"]`).
   Sim sample names must include a genome token (`mm10` or `hg38`), otherwise
   config generation and workflow initialization will fail fast.

5. Set the `CORE_NUM` variable to the number of CPU cores you want to use:

   CORE_NUM=64

   Replace `64` with the actual number of CPU cores available on your system.

6. Run the Snakemake workflow using the following command:

   snakemake --cores $CORE_NUM --profile snakemake_profile \
     --singularity-args "--bind /path/to/apabenchmark_final:/path/to/apabenchmark_final --cpus ${CORE_NUM}"

   This command will:
   - Use the specified number of CPU cores (`$CORE_NUM`).
   - Use the specified Snakemake profile (`snakemake_profile`).
   - Pass a project-root bind to Singularity/Apptainer.

7. Once the workflow finishes, the benchmarking results will be generated and saved in the designated output directory.

Note: Make sure you are in the directory where the benchmarking workflow is located before running the commands.
