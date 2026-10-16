version 1.0

task trimmomatic_task {
  meta {
    # task metadata
    description: "trims paired-end read data with trimmomatic tool"
  }
  input {
    # task inputs
    File read1
    File read2
    String docker = "quay.io/staphb/trimmomatic:0.39"
    # higher cpu and memory than fastq_scan
    Int cpu = 4
    Int memory = 8
    # trimmomatic specific parameters (different between exercise directions and solution)
    Int minlen = 75
    Int window_size = 10
    Int required_quality = 20
    Int threads = 4
  }
  command <<<
  # code block executed to run trimmomatic
  trimmomatic PE \
    -threads ~{threads} \
    ~{read1} ~{read2} \
    -baseout trimmed.fastq.gz \
    SLIDINGWINDOW:~{window_size}:~{required_quality} \
    MINLEN:~{minlen}
  >>>
  output {
    # task outputs
    File read1_trimmed = "trimmed_1P.fastq.gz"
    File read2_trimmed = "trimmed_2P.fastq.gz"
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
