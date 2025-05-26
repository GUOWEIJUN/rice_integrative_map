1.geneome preparation

bismark_genome_preparation --verbose genome

2.align

bismark --temp_dir ./1_align --multicore 16  --non_directional -o ./1_align/   --rg_id  $sample     --rg_sample $sample    --genome $genome  -1 $sample_R1.fq.gz         -2 $sample_R2.fq.gz

3.deduplicate remove

deduplicate_bismark --output_dir ./2_duplicate  ./bam/$sample.bam 

4.methylation_extractor

bismark_methylation_extractor --comprehensive --ignore_r2 2 -o  ./3_extractor  --gzip --multicore 16 --buffer_size 20G --bedGraph --cytosine_report      --genome_folder  $genome  $sample.deduplicated.bam  
