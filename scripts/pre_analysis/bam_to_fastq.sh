#!/bin/bash

for FILE in *.bam; do 
    mkdir s"$FILE" 
    samtools bam2fq $FILE > sample.fastq
    cat sample.fastq | grep '^@.*/1$' -A 3 --no-group-separator > SAMPLE_r1.fast
    cat sample.fastq | grep '^@.*/2$' -A 3 --no-group-separator > SAMPLE_r2.fast
    rm sample.fastq
    mv *.fastq s"$FILE" 
    kallisto quant -i Mus_musculusGRCh38.cdna.all.release-99.idx -t 20 -o s"$FIL
    cd s"$FILE"
    rm *.fastq
    cd ..
done

