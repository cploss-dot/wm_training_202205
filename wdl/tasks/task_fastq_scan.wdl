version 1.0

task fastq_scan_task {
  meta {
    # task metadata
    description: "Task to run fastq_scan with pairend read files"
  }
  input {
    # task inputs
    File read1
    File read2
    String docker = "staphb/fastq-scan:0.4.4"
    Int cpu = 2
    Int memory = 2
  }
  command <<<
    # code block executed 
    # ~curly_brackets for wdl variables 
    # same command as earlier part of exercise 
    # same result in new file 
    zcat ~{read1} | fastq-scan | jq .qc_stats.read_total > READ1_TOTAL_READS
    zcat ~{read2} | fastq-scan | jq .qc_stats.read_total > READ2_TOTAL_READS
  >>>
  output {
    # task outputs
    # takes files created in command block
    Int read1_total_reads = read_string("READ1_TOTAL_READS")
    Int read2_total_reads = read_string("READ2_TOTAL_READS")
  }
  runtime {
    # runtime environment
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
