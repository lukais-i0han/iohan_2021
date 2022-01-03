#/bin/bash

### Script to make the index used by Kallisto to do the pseudo-aligment. This index was used for 5XFAD and TauD35 pseudo-aligments

curl -O ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx Homo_sapiens.GRCh38.cdna.all.fa.gz

