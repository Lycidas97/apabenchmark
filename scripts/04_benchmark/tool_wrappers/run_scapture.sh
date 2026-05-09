# conda activate SCAPTURE_env
# export PATH="/path/to/apabenchmark_final/benchmark/bin/SCAPTURE-main":$PATH
export PATH="/scapture_tools/SCAPTURE-main":$PATH

# bamfile=`realpath $1`
# workdir=$2
# ncores=$4
# length=$5

while getopts b:o:g:l:t:c: flag
do
    case "${flag}" in
        b) bamfile=${OPTARG};;
        o) outputdir=${OPTARG};;
        g) genome=${OPTARG};;
        l) length=${OPTARG};;
        t) threads=${OPTARG};;
        c) cellbarcode=${OPTARG};;
    esac
done

if [ $genome == "mm10" ]; then
    genome_prefix="/scapture_annotation/mm10"
    genome_fa="/scapture_annotation/GRCm38.p6.genome.fa"
    polyaDB="/scapture_annotation/mm10.SupTab_KnownPASs_fourDBs.txt"
    species="mouse"

elif [ $genome == "hg38" ]; then
    genome_prefix="/scapture_annotation/hg38"
    genome_fa="/scapture_annotation/GRCh38.p14.genome.fa"
    polyaDB="/scapture_annotation/hg38.SupTab_KnownPASs_fourDBs.txt"
    species="human"
else
    echo "Error: genome not supported"
    exit 1
fi

if [ -z "$bamfile" ] || [ -z "$outputdir" ] || [ -z "$length" ] || [ -z "$genome" ]; then
    echo "Usage: $0 -b <bamfile> -o <outputdir> -g <genome> -l <length> -t <threads>"
    exit 1
fi

# for test
# genome_prefix="/path/to/apabenchmark_final/benchmark/env/scapture/mm10"
# genome_fa="/path/to/apabenchmark_final/benchmark/env/scapture/GRCm38.p6.genome.fa"
# polyaDB="/path/to/apabenchmark_final/benchmark/env/scapture/mm10.SupTab_KnownPASs_fourDBs.txt"
# species="mouse"

# bamfile="/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_gn1000_rep1.bam"
# cellbarcode="/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_gn1000_rep1.barcode_list.txt"
# outputdir="/path/to/apabenchmark_final/benchmark/result/scapture/test"
# length=120
# threads=8


if [ ! -d "$outputdir" ]; then
    mkdir -p $outputdir
fi

echo "Start"
scapture -m PAScall \
    -a $genome_prefix \
    -g $genome_fa \
    -b $bamfile \
    -l $length \
    -o $outputdir/tmp \
    -p $threads \
    --species $species \
    --polyaDB $polyaDB

# perl -alne '$,="\t";print @F[0..11] if $F[12] > 0 | $F[13] eq "positive";' .exonic.peaks.bed .intronic.peaks.evaluated.bed > pas.bed

perl -alne '$,="\t";print @F[0..11];' $outputdir/tmp.exonic.peaks.evaluated.bed $outputdir/tmp.intronic.peaks.evaluated.bed > $outputdir/tmp.PASquant.bed


scapture -m PASquant \
    -o $outputdir/tmp \
    -p $threads \
    -b $bamfile \
    --pas $outputdir/tmp.PASquant.bed \
    --celllist $cellbarcode

gzip -d $outputdir/tmp.KeepCell.UMIs.tsv.gz 

# awk 'BEGIN {FS="\t"; OFS="\t"} 
#     NR == 1 {
#         for (i = 2; i <= NF; i++) {
#             sub(/-1$/, "", $i)  # 移除列标题中末尾的-1
#             sub(/^X/, "", $i)  # 移除列标题中开头的X
#         }
#     }
#     {
#         print
#     }' $outputdir/tmp.KeepCell.UMIs.tsv | sed '1s/^[^\t]*\t//' > $outputdir/tmp.pas_counts_to_transpose.tsv

awk 'BEGIN {FS="\t"; OFS="\t"}
    NR == 1 {
        for (i = 2; i <= NF; i++) {
            gsub(/-1$/, "", $i)  # 移除列标题中末尾的-1
            len = length($i)
            count[len]++
            if (count[len] > maxCount) {
                maxCount = count[len]
                absLen = len
            }
        }
    }
    {
        print
    }' $outputdir/tmp.KeepCell.UMIs.tsv | sed '1s/^[^\t]*\t//' > $outputdir/tmp.pas_counts_to_transpose.tsv



python -c "import sys, pandas as pd; \
df = pd.read_csv(sys.stdin, sep='\t'); \
df = df.T; \
mode_length = pd.Series([len(str(x)) for x in df.index]).mode()[0]; \
df.index = [r[1:] if len(r)>mode_length and r.startswith('X') else r for r in df.index]; \
df.to_csv(sys.stdout, sep='\t')" < $outputdir/tmp.pas_counts_to_transpose.tsv | sed '1s/^[^\t]*\t//' > $outputdir/pas_counts.tsv

awk 'BEGIN { OFS="\t" }
{
    if ($6 == "+") {
        print $1, $3, $3 + 1, $4, $5, $6
    } else if ($6 == "-") {
        print $1, $2 - 1, $2, $4, $5, $6
    } else {
        print $0
    }
}' $outputdir/tmp.exonic.peaks.bed > $outputdir/tmp.unsorted.pas.bed
awk 'BEGIN { OFS="\t" }
{
    if ($6 == "+") {
        print $1, $3, $3 + 1, $4, $5, $6
    } else if ($6 == "-") {
        print $1, $2 - 1, $2, $4, $5, $6
    } else {
        print $0
    }
}' $outputdir/tmp.intronic.peaks.bed >> $outputdir/tmp.unsorted.pas.bed

bedtools sort -i $outputdir/tmp.unsorted.pas.bed | awk 'BEGIN{FS=OFS="\t"} {sub(/\|.*/, "", $4); if ($4 in count) {count[$4]++} else {count[$4]=1}; $4=$4"-"count[$4]} 1' > $outputdir/pas.bed

rm -rf $outputdir/_tmp
rm -rf $outputdir/tmp*

# exec /bin/bash -c 'source /scapture/bin/activate &&         /run_scapture.sh -b /home/singularity_project/apabenchmark/benchmark/../sim_bam/Visium_mouse_brain_V19L01-041-C1_pas1_gn5000_rep2.bam -o /home/singularity_project/apabenchmark/benchmark/sim_bam_result/scapture/Visium_mouse_brain_V19L01-041-C1_pas1_gn5000_rep2 -g mm10 -l 120 -t 4 > /home/singularity_project/apabenchmark/benchmark/sim_bam_result/scapture/Visium_mouse_brain_V19L01-041-C1_pas1_gn5000_rep2/log.txt'