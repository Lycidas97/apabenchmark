config_file="./snakemake_profile/sample.yaml"
rawbam_dir="../../data/int_data/bam_to_detect_pas"
echo "raw_sample:" > $config_file
for file in $rawbam_dir/*.bam.bai; do
    bamfile=${file%.bai}
    filename=$(basename "$file")
    filename=${filename%.bai}
    sample=${filename%.bam}
    read_length=$(samtools view $bamfile | head -10000 | awk '{print length($10)}' | sort | uniq -c | sort -nrk1,1 | awk 'NR==1 {print $2}')
    echo "  $sample:" >> $config_file
    echo "    read_length: $read_length" >> $config_file
done

simbam_dir="../../data/sim_data/bam"
echo "sim_sample:" >> $config_file
for file in $simbam_dir/*.bam.bai; do
    bamfile=${file%.bai}
    filename=$(basename "$file")
    filename=${filename%.bai}
    sample=${filename%.bam}
    read_length=$(samtools view $bamfile | head -10000 | awk '{print length($10)}' | sort | uniq -c | sort -nrk1,1 | awk 'NR==1 {print $2}')
    echo "  $sample:" >> $config_file
    echo "    read_length: $read_length" >> $config_file
done

pf_bam_dir="../../data/sim_data/bam_for_cprsb"
echo "pf_sample:" >> $config_file
for file in $pf_bam_dir/*.bam.bai; do
    bamfile=${file%.bai}
    filename=$(basename "$file")
    filename=${filename%.bai}
    sample=${filename%.bam}
    read_length=$(samtools view $bamfile | head -10000 | awk '{print length($10)}' | sort | uniq -c | sort -nrk1,1 | awk 'NR==1 {print $2}')
    echo "  $sample:" >> $config_file
    echo "    read_length: $read_length" >> $config_file
done