#!/bin/bash -ue
#==========================================================================================
# Assemble based on runtype
#==========================================================================================
mkdir supplemental
if [ "false" == "true" ]; then
    # Skip assembly and use provided assembly
    echo "Using provided assembly for SRR2838702 without re-assembly." > supplemental/assembly-info.txt
    gzip -cd null > SRR2838702.fna

    # passthrough correct reads for downstream analysis
    if [[ "paired-end" == "short_polish" ]] || [[ "paired-end" == "ont" ]]; then
        cp -L null SRR2838702_ONT.fastq.gz
        touch supplemental/ont.txt
    elif [[ "false" == "true" ]]; then
        cp -L null SRR2838702_SE.fastq.gz
        touch supplemental/illumina.txt
    else
        cp -L staging/reads/r1/SRR2838702_R1.fastq.gz SRR2838702_R1.fastq.gz
        cp -L staging/reads/r2/SRR2838702_R2.fastq.gz SRR2838702_R2.fastq.gz
        touch supplemental/illumina.txt
    fi
elif [[ "paired-end" == "hybrid" || "false" == "true" ]]; then
    #======================================================================================
    # Unicycler Assembler
    #
    # Unicycler is used for hybrid assemblies (short + long reads) or when explicitly
    # requested via the --use_unicycler parameter.
    #======================================================================================
    unicycler \
        -1 staging/reads/r1/SRR2838702_R1.fastq.gz -2 staging/reads/r2/SRR2838702_R2.fastq.gz \
        --keep 1 --min_fasta_length 500 --mode normal --min_component_size 1000 --min_dead_end_size 1000 \
         \
        -o supplemental/ \
        --threads 4
    sed -r 's/^>([0-9]+)(.*)/>SRR2838702_\1\2/' supplemental/assembly.fasta > SRR2838702.fna
    mv supplemental/assembly.fasta supplemental/unicycler-unpolished.fasta
    mv supplemental/assembly.gfa supplemental/unicycler-unpolished.gfa

    # Illumina FASTQs should be used downstream
    cp -L staging/reads/r1/SRR2838702_R1.fastq.gz SRR2838702_R1.fastq.gz
    cp -L staging/reads/r2/SRR2838702_R2.fastq.gz SRR2838702_R2.fastq.gz
    touch supplemental/illumina.txt
elif [[ "paired-end" == "ont" || "paired-end" == "short_polish" ]]; then
    #======================================================================================
    # Dragonflye Assembler
    #
    # Dragonflye is used for long read assemblies (ONT) with optional short read polishing.
    #======================================================================================
    if ! dragonflye \
        --reads null \
        --gsize 358242 \
        --outdir supplemental \
        --assembler flye --polypolish 1 --minlen 500 --mincov 2 --force --keepfiles --depth 0 --minreadlen 0 --minquality 0 --racon 1 --medaka 0 \
        --namefmt "SRR2838702_%05d" \
        --cpus 4 \
        --ram 7 \
        --noreorient; then

        # Check if error is due to no contigs
        if grep "has zero contigs" supplemental/dragonflye.log; then
            touch supplemental/contigs.fa
        else
            exit 1
        fi
    fi
    mv supplemental/contigs.fa SRR2838702.fna

    # ONT FASTQ should be used downstream
    cp -L null SRR2838702_ONT.fastq.gz
    touch supplemental/ont.txt
else
    #======================================================================================
    # Shovill Assembler
    #
    # Shovill is used for short read assemblies (Illumina PE or SE).
    #======================================================================================
    if ! shovill --R1 staging/reads/r1/SRR2838702_R1.fastq.gz --R2 staging/reads/r2/SRR2838702_R2.fastq.gz \
        --gsize 358242 \
        --outdir supplemental \
        --assembler skesa --minlen 500 --mincov 2 --force --keepfiles --depth 0 --noreadcorr \
        --namefmt "SRR2838702_%05d" \
        --cpus 4 \
        --ram 7; then

        # Check if error is due to no contigs
        if grep "has zero contigs" supplemental/shovill.log; then
            touch supplemental/contigs.fa
        else
            exit 1
        fi
    fi
    mv supplemental/contigs.fa SRR2838702.fna

    # Rename Graphs
    if [ -f "supplemental/contigs.gfa" ]; then
        mv supplemental/contigs.gfa supplemental/skesa-unpolished.gfa
    elif [ -f "supplemental/contigs.fastg" ]; then
        mv supplemental/contigs.fastg supplemental/skesa-unpolished.gfa
    elif [ -f "supplemental/contigs.LastGraph" ]; then
        mv supplemental/contigs.LastGraph supplemental/skesa-unpolished.gfa
    fi

    if [ -f "supplemental/flye-info.txt" ]; then
        mv supplemental/flye-info.txt supplemental/flye.log
    fi

    # Illumina FASTQs should be used downstream
    if [ "false" == "false" ]; then
        cp -L staging/reads/r1/SRR2838702_R1.fastq.gz SRR2838702_R1.fastq.gz
        cp -L staging/reads/r2/SRR2838702_R2.fastq.gz SRR2838702_R2.fastq.gz
    else
        cp -L null SRR2838702_SE.fastq.gz
    fi
    touch supplemental/illumina.txt
fi

#==========================================================================================
# Assembly Quality Check
#==========================================================================================
TOTAL_CONTIGS=`grep -c "^>" SRR2838702.fna || true`
if [ "${TOTAL_CONTIGS}" -gt 0 ]; then
    assembly-scan SRR2838702.fna --prefix SRR2838702 > SRR2838702.tsv
    TOTAL_CONTIG_SIZE=$(cut -f 3 SRR2838702.tsv | tail -n 1)
    if [ "${TOTAL_CONTIG_SIZE}" -lt null ]; then
        mv SRR2838702.fna SRR2838702-error.fna
        mv SRR2838702.tsv SRR2838702-error.tsv
        echo "SRR2838702 assembled size [${TOTAL_CONTIG_SIZE} bp] is less than the minimum allowed genome
                size [null bp]. If this is unexpected, please investigate SRR2838702 to
                determine a cause [e.g. metagenomic, contaminants, etc...] for the poor assembly.
                Otherwise, adjust the --min_genome_size parameter to fit your need. Further assembly
                based analysis of SRR2838702 will be discontinued." |             sed 's/^\s*//' > SRR2838702-assembly-error.txt
    elif [ "${TOTAL_CONTIG_SIZE}" -gt null ]; then
        mv SRR2838702.fna SRR2838702-error.fna
        mv SRR2838702.tsv SRR2838702-error.tsv
        echo "SRR2838702 assembled size [${TOTAL_CONTIG_SIZE} bp] is greater than the maximum allowed genome
                size [null bp]. If this is unexpected, please investigate SRR2838702 to
                determine a cause [e.g. metagenomic, contaminants, etc...] for the poor assembly.
                Otherwise, adjust the --max_genome_size parameter to fit your need. Further assembly
                based analysis of SRR2838702 will be discontinued." |             sed 's/^\s*//' > SRR2838702-assembly-error.txt
    fi
else
    mv SRR2838702.fna SRR2838702-error.fna
    echo "SRR2838702 assembled successfully, but 0 contigs were formed. Please investigate
            SRR2838702 to determine a cause [e.g. metagenomic, contaminants, etc...] for this
            outcome. Further assembly-based analysis of SRR2838702 will be discontinued." |         sed 's/^\s*//' > SRR2838702-assembly-error.txt
fi

#==========================================================================================
# Cleanup
#==========================================================================================
if [ "false" == "false" ]; then
    # Remove intermediate files
    rm -rfv supplemental/shovill.bam* \
            supplemental/shovill-se.bam* \
            supplemental/flash.extendedFrags* \
            supplemental/flash.notCombined* \
            supplemental/skesa.fasta* \
            supplemental/*.fq.gz \
            supplemental/00*.gfa \
            supplemental/pilon_polish* \
            supplemental/flye/ \
            supplemental/flye.fasta* \
            supplemental/raven/ \
            supplemental/raven.fasta* \
            supplemental/raven.cereal \
            supplemental/miniasm/ \
            supplemental/miniasm.fasta* \
            supplemental/spades/ \
            supplemental/spades.fasta* \
            supplemental/megahit/ \
            supplemental/megahit.fasta* \
            supplemental/velvet.fasta* \
            supplemental/velvet/ \
            supplemental/assembly.fasta
fi

#==========================================================================================
# Compress and move final outputs
#==========================================================================================
if [[ "false" == "false" ]]; then
    # Compress based on matched extensions
    pigz -n --best -p 4 SRR2838702*.fna
    find supplemental/ -type f |             grep -E "\.fna$|\.fasta$|\.fa$|\.gfa$" |             xargs -I {} pigz -n --best -p 4 {}
fi
find supplemental/ -maxdepth 1 -name "*.log" | xargs -I {} mv {} ./

# Capture versions (common tools available on all platforms)
cat <<-END_VERSIONS > versions.yml
"ASSEMBLER:ASSEMBLER_MODULE":
    assembly-scan: $(echo $(assembly-scan --version 2>&1) | sed 's/assembly-scan //')
    bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
    flash: $(echo $(flash --version 2>&1) | sed 's/^.*FLASH v//;s/ .*$//')
    megahit: $(echo $(megahit --version 2>&1) | sed 's/MEGAHIT v//')
    miniasm: $(echo $(miniasm -V))
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    pilon: $(echo $(pilon --version 2>&1) | sed 's/^.*Pilon version //;s/ .*$//')
    racon: $(echo $(racon --version 2>&1) | sed 's/v//')
    samclip: $(echo $(samclip --version 2>&1) | sed 's/^.*samclip //')
    samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
    shovill: $(echo $(shovill --version 2>&1) | sed 's/^.*shovill //')
    shovill-se: $(echo $(shovill-se --version 2>&1) | sed 's/^.*shovill-se //')
    skesa: $(echo $(skesa --version 2>&1) | sed 's/^.*SKESA //;s/ .*$//')
    spades.py: $(echo $(spades.py --version 2>&1) | sed 's/SPAdes genome assembler v//')
    velvetg: $(echo $(velvetg 2>&1) | sed 's/^.*Version //;s/ .*$//')
    velveth: $(echo $(velveth 2>&1) | sed 's/^.*Version //;s/ .*$//')
    unicycler: $(echo $(unicycler --version 2>&1) | sed 's/^.*Unicycler v//;s/ .*$//')
END_VERSIONS

# Linux-only tools (long-read assemblers not available on macOS)
if [[ "$OSTYPE" != "darwin"* ]]; then
cat <<-END_LINUX >> versions.yml
    any2fasta: $(echo $(any2fasta -v 2>&1) | sed 's/^.*any2fasta //')
    dragonflye: $(echo $(dragonflye --version 2>&1) | sed 's/^.*dragonflye //' )
    flye: $(echo $(flye --version))
    medaka: $(echo $(medaka --version 2>&1) | sed 's/medaka //')
    minimap2: $(echo $(minimap2 --version))
    nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
    rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
    raven: $(echo $(raven --version))
END_LINUX
fi
