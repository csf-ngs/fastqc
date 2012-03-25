###private fork of fastqc (0.10.0)

added a read argument for unaligned bam files.
probably breaks other file types, but I don't care.

#####call with:
   fastqc-csf/fastqc --read {d} file.bam

   where {d} is:
   
   * 0: unaligned bam file with one read only (4 in flags)
   * 1: unaligned bam file with paired end reads. Takes only 1. read
   * 2: unaligned bam file with paired end reads. Takes only 2. read

#####build:
   1. ant deploy.
   2. copy the fastqc-csf folder to your destination.



