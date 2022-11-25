## cram_to_fastq

### dependencies
samtools/1.9 
biobambam/2.0.87--1 
nextflow/22.04.5-5708

### usage
nextflow run cram_to_fastq.nf --manifest <manifest> --reference <reference> -profile sanger_lsf

### notes
Designed to run on the sanger farm5 HPC system
There is an example manifest in this repo
