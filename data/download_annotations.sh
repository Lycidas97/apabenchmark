#!/bin/bash
mkdir -p raw_data/annotations

wget -O "raw_data/annotations/atlas.clusters.2.0.GRCm38.96.bed" "https://bis.zju.edu.cn/nextcloud/s/apabenchmark_annotation/download?path=%2F&files=atlas.clusters.2.0.GRCm38.96.bed"
wget -O "raw_data/annotations/gencode.vM25.polyAs.gtf" "https://bis.zju.edu.cn/nextcloud/s/apabenchmark_annotation/download?path=%2F&files=gencode.vM25.polyAs.gtf"
wget -O "raw_data/annotations/mm10.3utr.sorted.bed" "https://bis.zju.edu.cn/nextcloud/s/apabenchmark_annotation/download?path=%2F&files=mm10.3utr.sorted.bed"
wget -O "raw_data/annotations/gencode.vM25.annotation.bed" "https://bis.zju.edu.cn/nextcloud/s/apabenchmark_annotation/download?path=%2F&files=gencode.vM25.annotation.bed"
wget -O "raw_data/annotations/mouse.PAS.txt" "https://bis.zju.edu.cn/nextcloud/s/apabenchmark_annotation/download?path=%2F&files=mouse.PAS.txt"
