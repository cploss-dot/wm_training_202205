#!/bin/bash -ue
#==========================================================================================
# Helper functions
#===========================================================================================
check_gzip_or_empty() {
    local fq="$1"
    if ! gzip -t "${fq}"; then
        # Failed Gzip test
        echo "SRR2838702 FASTQs failed Gzip tests. Please check the input FASTQs. Further analysis is discontinued." > "SRR2838702-gzip-error.txt"
        return 1
    fi
    # Empty FASTQ file
    echo "SRR2838702 FASTQs are empty. Please check the input FASTQs. Further analysis is discontinued." > "SRR2838702-empty-error.txt"
    return 1
}

MERGED="multiple-read-sets-merged.txt"
merge_runs() {
    local label="$1"
    local pattern="$2"
    local output="$3"

    echo "${label}:" >> ${MERGED}
    find -name "${pattern}" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"	"$9}' >> ${MERGED}
    find -name "${pattern}" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > "${output}"
    echo "Merged ${label}:" >> ${MERGED}
    ls -l "${output}" | awk '{print $5"	"$9}' >> ${MERGED}
}

write_download_error() {
    local source="$1"  # "SRA and ENA" or "NCBI Assembly"
    local error_file="$2"
    echo "Unable to download SRR2838702 from ${source} 3 times. This may or may not be a temporary connection issue. Rather than stop the whole Bactopia run, further analysis of SRR2838702 will be discontinued." > "${error_file}"
}

validate_single() {
    local fq="$1"
    gzip -cd "${fq}" | fastq-scan > r1.json

    if [[ -s r1.json ]]; then
        bactopia-check-fastqs --fq1 r1.json ${OPTS} || { rm -f r1.json; return 1; }
    else
        check_gzip_or_empty "${fq}"
        rm -f r1.json
        return 1
    fi

    rm -f r1.json
    return 0
}

#==========================================================================================
# Gather or simulate reads based on runtype
#==========================================================================================
mkdir -p fastqs assembly
case "paired-end" in
    paired-end)
        # Paired-End Reads
        cp -L staging/r1/001-r1 fastqs/SRR2838702_R1.fastq.gz
        cp -L staging/r2/001-r2 fastqs/SRR2838702_R2.fastq.gz
        ;;
    single-end)
        # Single-End Reads
        cp -L null fastqs/SRR2838702_SE.fastq.gz
        ;;
    ont)
        # Long reads (Nanopore)
        cp -L null fastqs/SRR2838702_ONT.fastq.gz
        ;;
    hybrid|short_polish)
        # Paired-End Reads + Long Reads
        cp -L staging/r1/001-r1 fastqs/SRR2838702_R1.fastq.gz
        cp -L staging/r2/001-r2 fastqs/SRR2838702_R2.fastq.gz
        cp -L null fastqs/SRR2838702_ONT.fastq.gz
        ;;
    merge-pe|hybrid-merge-pe|short_polish-merge-pe)
        # Merged Paired-End Reads
        echo "This sample had reads merged." > ${MERGED}
        merge_runs "R1" "*-r1" "fastqs/SRR2838702_R1.fastq.gz"
        merge_runs "R2" "*-r2" "fastqs/SRR2838702_R2.fastq.gz"
        if [ "paired-end" == "hybrid-merge-pe" ] || [ "paired-end" == "short_polish-merge-pe" ]; then
            cp -L null fastqs/SRR2838702_ONT.fastq.gz
        fi
        ;;
    merge-se)
        # Merged Single-End Reads
        echo "This sample had reads merged." > ${MERGED}
        merge_runs "SE" "*-se" "fastqs/SRR2838702_SE.fastq.gz"
        ;;
    sra_accession|sra_accession_ont)
        if [ "1" == "3" ]; then
            write_download_error "SRA and ENA" "SRR2838702-fastq-download-error.txt"
            exit
        else
            # Download accession from ENA/SRA
            fastq-dl \
                --accession SRR2838702 \
                --provider SRA \
                --cpus 4 \
                --outdir fastqs/ \
                --prefix SRR2838702 \
                --group-by-experiment

            # Rename single FASTQ files based on runtype
            if [ -f "fastqs/SRR2838702.fastq.gz" ]; then
                if [ "paired-end" == "sra_accession_ont" ]; then
                    # ONT reads
                    mv fastqs/SRR2838702.fastq.gz fastqs/SRR2838702_ONT.fastq.gz
                else
                    # Single-End
                    mv fastqs/SRR2838702.fastq.gz fastqs/SRR2838702_SE.fastq.gz
                fi
            fi
        fi
        ;;
    assembly_accession|assembly)
        if [ "paired-end" == "assembly_accession" ]; then
            if [ "1" == "3" ]; then
                write_download_error "NCBI Assembly" "SRR2838702-assembly-download-error.txt"
                exit
            fi

            # Verify Assembly accession
            bactopia-check-assembly-accession SRR2838702 > accession.txt 2> check-assembly-accession.txt

            if [ ! -s "accession.txt" ]; then
                mv check-assembly-accession.txt SRR2838702-assembly-accession-error.txt
                exit
            else
                rm check-assembly-accession.txt
            fi

            # Download from NCBI Datasets and simulate reads
            genome-dl --formats fasta --outdir ./ --cpus 4 --max-attempts 3 \
                      --prefix SRR2838702 --accession SRR2838702

            # genome-dl names the assembly by its resolved accession. The assembly version
            # (e.g., GCF_000005845.2 --> .2) is removed for consistency.
            find . -name "*SRR2838702*.fna.gz" | xargs -I {} mv {} assembly/
            rename 's/(GC[AF]_\d+).*/$1.fna.gz/' assembly/*
            gzip -cd assembly/SRR2838702.fna.gz > SRR2838702-art.fna
        elif [ "paired-end" == "assembly" ]; then
            if [ "true" == "true" ]; then
                gzip -cd null > SRR2838702-art.fna
            else
                cat null > SRR2838702-art.fna
            fi
            # Keep the original assembly
            cp -L null assembly/SRR2838702.fna.gz
        fi

        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov 150 -ir 0 -ir2 0 -dr 0 -dr2 0 -rs null\
                     -na -qL 33 -qU 40 -o SRR2838702_R --id SRR2838702 -i SRR2838702-art.fna

        mv SRR2838702_R1.fq fastqs/SRR2838702_R1.fastq
        mv SRR2838702_R2.fq fastqs/SRR2838702_R2.fastq
        pigz -p 4 --fast fastqs/*.fastq
        rm SRR2838702-art.fna
        ;;
    *)
        # This should not happen, in case it does fail the whole pipeline
        echo "Runtype paired-end is not recognized. Further analysis of SRR2838702 will be discontinued." > SRR2838702-runtype-error.txt
        exit 1
        ;;
esac

#==========================================================================================
# Validate input FASTQs
#==========================================================================================
if [ "false" == "false" ]; then
    ERROR=0
    # Check paired-end reads have same read counts
    OPTS="--sample SRR2838702 --min_basepairs 2241820 --min_reads 7472 --min_proportion 0.5 --runtype paired-end"
    if [ -f  "fastqs/SRR2838702_R2.fastq.gz" ]; then
        # Paired-end
        gzip -cd fastqs/SRR2838702_R1.fastq.gz | fastq-scan > r1.json
        gzip -cd fastqs/SRR2838702_R2.fastq.gz | fastq-scan > r2.json
        if ! reformat.sh in1=fastqs/SRR2838702_R1.fastq.gz in2=fastqs/SRR2838702_R2.fastq.gz qin=auto out=/dev/null 2> SRR2838702-paired-end-error.txt; then
            echo "SRR2838702 FASTQs contains an error. Please check the input FASTQs. Further analysis is discontinued." >> SRR2838702-paired-end-error.txt
            ERROR=1
        else
            rm -f SRR2838702-paired-end-error.txt
        fi

        if [[ -s r1.json ]] && [[ -s r2.json ]]; then
            bactopia-check-fastqs --fq1 r1.json --fq2 r2.json ${OPTS} || ERROR=1
        else
            check_gzip_or_empty "fastqs/SRR2838702_R1.fastq.gz" || ERROR=1
            check_gzip_or_empty "fastqs/SRR2838702_R2.fastq.gz" || ERROR=1
        fi
        rm -f r1.json r2.json
    elif [ -f "fastqs/SRR2838702_SE.fastq.gz" ]; then
        # Single-end
        validate_single "fastqs/SRR2838702_SE.fastq.gz" || ERROR=1
    elif [ -f "fastqs/SRR2838702_ONT.fastq.gz" ]; then
        # Long reads
        validate_single "fastqs/SRR2838702_ONT.fastq.gz" || ERROR=1
    fi

    # Failed validations so, let's keep them from continuing
    if [ "${ERROR}" -eq "1" ]; then
        mv fastqs/ failed-tests-fastqs/
    fi
fi

#==========================================================================================
# Determine the paired/single-end status
#==========================================================================================
IS_PAIRED="false"
if [ -f "fastqs/SRR2838702_R2.fastq.gz" ] && [ "paired-end" != "short_polish" ]; then
    IS_PAIRED="true"
fi

#==========================================================================================
# Create metadata file
#==========================================================================================
printf "sample\truntype\toriginal_runtype\tis_paired\tis_compressed\tspecies\tgenome_size\n" > SRR2838702-meta.tsv
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "SRR2838702" "paired-end" "paired-end" "$IS_PAIRED" "true" "null" "358242" >> SRR2838702-meta.tsv

# Cleanup

cat <<-END_VERSIONS > versions.yml
"GATHER:GATHER_MODULE":
    art: $(echo $(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*$//')
    fastq-dl: $(echo $(fastq-dl --version 2>&1) | sed 's/fastq-dl, version //')
    fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
    genomedl: $(echo $(genome-dl --version 2>&1) | sed 's/.*version //')
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
END_VERSIONS
