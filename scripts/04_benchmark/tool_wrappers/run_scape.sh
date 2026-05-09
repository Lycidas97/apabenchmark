#!/bin/bash

# 处理输入参数
scape_py_path="/scape/SCAPE"
# scape_py_path="/path/to/apabenchmark_final/benchmark/bin/SCAPE"
while getopts b:r:o:t:c:g: flag
do
    case "${flag}" in
        b) bamfile=${OPTARG};;
        r) bedfile=${OPTARG};;
        o) outputdir=${OPTARG};;
        t) threads=${OPTARG};;
        c) cellbarcode=${OPTARG};;
        g) tag=${OPTARG};;

    esac
done

# exit if bam bed outputdir threads not provided
if [ -z "$bamfile" ] || [ -z "$bedfile" ] || [ -z "$outputdir" ] || [ -z "$tag" ]; then
    echo "Usage: $0 -b <bamfile> -r <bedannotation> -o <outputdir> -t <threads> -g <tags>"
    exit 1
fi

# if use relative path, transform to absolute path
if [[ $bamfile != /* ]]; then
    bamfile=$PWD/$bamfile
fi
if [[ $outputdir != /* ]]; then
    outputdir=$PWD/$outputdir
fi
if [ -z "$threads" ]; then
    threads=1
fi


mkdir -p $outputdir

cd $scape_py_path

# generate cellbarcode if not provided
if [ -z "$cellbarcode" ]; then
    cellbarcode=$outputdir/cellbarcode.tsv
    samtools view $bamfile | awk -F '\t' '{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) print substr($i,6)}' | sort | uniq > $cellbarcode
fi

# if [[ $cellbarcode != /* ]]; then
#     cellbarcode=$PWD/$cellbarcode
# fi

# run scape
echo "Start Running SCAPE, Check Progress in $outputdir/scape.log"
python main.py apamix \
--bed $bedfile \
--bam $bamfile \
--out $outputdir \
--cores $threads \
--cb $cellbarcode \
--tag $tag &> $outputdir/scape.log

echo "Finish Running SCAPE, parsing output"

# wait until the output file is generated
max_try=10
while [ ! -f $outputdir/pasite.csv.gz ] && [ $max_try -gt 0 ]; do
    sleep 1
    max_try=$((max_try-1))
done

if [ ! -f $outputdir/pasite.csv.gz ]; then
    touch $outputdir/pas.bed
    touch $outputdir/pas_counts.tsv
else
    gzip -d $outputdir/pasite.csv.gz 
    python -c "import sys, pandas as pd; df = pd.read_csv(sys.stdin, index_col=0); df.T.to_csv(sys.stdout, sep='\t')" < $outputdir/pasite.csv | sed '1s/^[^\t]*\t//' > $outputdir/pas_counts.tsv

    awk -F, 'NR > 1 {
        split($1, a, "[:]"); 
        strand = a[4];
        position = a[2];
        if (strand == "+") {
            start = position;
            end = position + 1;
        } else if (strand == "-") {
            start = position - 1;
            end = position;
        }
        print a[1]"\t"start"\t"end"\t"$1"\t"position"\t"strand
    }' $outputdir/pasite.csv > $outputdir/pas.bed

    find $outputdir -mindepth 1 -type d -exec rm -rf {} +
fi

echo "Done"