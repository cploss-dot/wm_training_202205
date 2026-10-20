#!/bin/bash -ue
#==========================================================================================
# Helper Functions
#==========================================================================================
is_assembly()        { [[ "false" == "true" ]]; }
is_hybrid()          { [[ "paired-end" =~ ^(hybrid|short_polish)$ ]]; }
is_ont()             { [[ "paired-end" == "ont" ]]; }
is_paired()          { [[ "false" == "false" ]]; }
is_illumina_se()     { [[ "false" == "true" ]] && ! is_ont; }
has_ont_reads()      { is_ont || is_hybrid; }
has_illumina_reads() { is_paired || is_illumina_se; }
fastq_scan()         { gzip -cd "$1" | fastq-scan -g 358242 > "$2"; }
compress()           { pigz -p 4 -c -n "$1" > "$2"; }

run_nanoplot() {
    local input="$1"
    local suffix="$2"
    local outdir="supplemental/SRR2838702-${suffix}"

    mkdir -p "${outdir}"
    NanoPlot              --threads 4             --fastq "${input}"             --outdir "${outdir}/"             --prefix "SRR2838702-${suffix}_"

    cp "${outdir}/SRR2838702-${suffix}_NanoPlot-report.html" "supplemental/SRR2838702-${suffix}_NanoPlot-report.html"
    tar -cvf - "${outdir}/" | pigz --best -p 4 > "supplemental/SRR2838702-${suffix}_NanoPlot.tar.gz"
    rm -rf "${outdir}"
}

report_empty_error() {
    echo "After ${1}, SRR2838702 FASTQs are empty. Please check the input FASTQs.
        Further analysis is discontinued." | sed 's/^ *//' >> "SRR2838702-${2}-qc-error.txt"
    ERROR=1
}

#==========================================================================================
# QC Initialization
#==========================================================================================
mkdir assembly supplemental
ERROR=0

if is_assembly; then
    # Copy the assembly over to assembly
    cp null assembly/SRR2838702.fna.gz
fi

#==========================================================================================
# Begin QC
#==========================================================================================
if [[ "false" == "true" ]]; then
    #======================================================================================
    # User skipped QC, just copy input reads to output
    #======================================================================================
    echo "Sequence QC was skipped for SRR2838702" > SRR2838702-qc-skipped.txt
    if is_paired; then
        cp staging/r1/SRR2838702_R1.fastq.gz SRR2838702_R1.fastq.gz
        cp staging/r2/SRR2838702_R2.fastq.gz SRR2838702_R2.fastq.gz
    fi

    if has_ont_reads; then
        cp null SRR2838702_ONT.fastq.gz
    fi

    if is_illumina_se; then
        cp staging/r1/SRR2838702_R1.fastq.gz SRR2838702_SE.fastq.gz
    fi
else
    #======================================================================================
    # Nanopore Read QC
    #======================================================================================
    if has_ont_reads; then
        NANOQ_FQ="null"
        if [[ "false" == "true" ]]; then
            # Remove Adapters
            NANOQ_FQ="adapter-ont.fq"
            porechop                     --input null                      --format fastq                     --threads 4 > adapter-ont.fq
        fi

        # Quality filter
        nanoq                 --min-len 1000                 --min-qual 0                 --input ${NANOQ_FQ} 1> filt-ont.fq
    fi

    #======================================================================================
    # Illumina Read QC
    #
    # Defaults to using Fastp unless BBMap is specified (--use_bbmap)
    #
    # ***** Note: Any errors in this block will set ERROR=1 and skip remaining steps *****
    #======================================================================================
    if has_illumina_reads; then
        if [[ "false" == "true" ]]; then
            # Use BBMap for cleaning reads
            if is_paired; then
                # Make sure paired-end reads have matching IDs
                repair.sh                         in=staging/r1/SRR2838702_R1.fastq.gz                         in2=staging/r2/SRR2838702_R2.fastq.gz                         out=repair-r1.fq                         out2=repair-r2.fq                         outs=repair-singles.fq                         ain=f

                if [ ! -s repair-r1.fq ]; then
                    report_empty_error "validating read pairs" "paired-match"
                fi
            else
                # No need to validate single-end read IDs
                gunzip -c staging/r1/SRR2838702_R1.fastq.gz > repair-r1.fq
            fi

            if [ "${ERROR}" -eq "0" ]; then
                # Remove Adapters
                bbduk.sh -Xmx8160437862                         in=repair-r1.fq out=adapter-r1.fq in2=repair-r2.fq out2=adapter-r2.fq                         ref=adapters                         k=23                         ktrim=r                         mink=11                         hdist=1                         tpe=t                         tbo=t                         threads=4                         ftm=5                         qin=auto ordered=t 

                if [ ! -s adapter-r1.fq ]; then
                    report_empty_error "adapter removal" "adapter"
                fi
            fi

            if [ "${ERROR}" -eq "0" ]; then
                # Remove PhiX
                bbduk.sh -Xmx8160437862                         in=adapter-r1.fq out=phix-r1.fq in2=adapter-r2.fq out2=phix-r2.fq                         ref=phix                         k=31                         hdist=1                         tpe=t                         tbo=t                         qtrim=rl                         trimq=6                         minlength=35                         minavgquality=10                         qin=auto qout=33                         tossjunk=t                         threads=4                         ordered=t 

                if [ ! -s phix-r1.fq ]; then
                    report_empty_error "PhiX removal" "phix"
                fi
            fi
        else
            # QC with fastp
            fastp                     --in1 staging/r1/SRR2838702_R1.fastq.gz --out1 filt-r1.fq --in2 staging/r2/SRR2838702_R2.fastq.gz --out2 filt-r2.fq --detect_adapter_for_pe                     --thread 4                     --json supplemental/SRR2838702.fastp.json                     --html supplemental/SRR2838702.fastp.html  2> SRR2838702-fastp.log
        fi
    fi

    # Error Correction
    if [ "${ERROR}" -eq "0" ]; then
        if [[ "false" == "true" ]]; then
            if [ "false" == "false" ] && [ "358242" -gt "0" ]; then
                lighter -od . -r phix-r1.fq -r phix-r2.fq -K 31 358242 -maxcor 1 -zlib 0 -t 4
                mv phix-r1.cor.fq filt-r1.fq
                if is_paired; then
                    mv phix-r2.cor.fq filt-r2.fq
                fi
            else
                echo "Skipping error correction"
                ln -s phix-r1.fq filt-r1.fq
                if is_paired; then
                    ln -s phix-r2.fq filt-r2.fq
                fi
            fi
        elif [[ "paired-end" == "ont" ]]; then
            echo "Skipping error correction. Have a recommended ONT error corrector? Let me know!"
        else
            echo "Skipping error correction since fastp was used"
        fi
    fi

    # Reduce Coverage
    if [ "${ERROR}" -eq "0" ]; then
        TOTAL_BP=$(( 100*358242 ))
        if (( ${TOTAL_BP} > 0 )); then
            if has_ont_reads; then
                rasusa reads                         -c 100                         -g 358242                         -s 42                         filt-ont.fq 1> subsample-ont.fq
            fi

            if has_illumina_reads; then
                if [ -f filt-r1.fq ]; then
                    rasusa reads                             -c 100                             -g 358242                             -s 42                             -o subsample-r1.fq -o subsample-r2.fq filt-r1.fq filt-r2.fq
                fi
            fi
        else
            echo "Skipping coverage reduction"
            if has_ont_reads; then
                ln -s filt-ont.fq subsample-ont.fq
            fi

            if has_illumina_reads; then
                ln -s filt-r1.fq subsample-r1.fq
                if is_paired; then
                    ln -s filt-r2.fq subsample-r2.fq
                fi
            fi
        fi
    fi

    # Compress
    if [ "${ERROR}" -eq "0" ]; then
        if is_paired; then
            compress subsample-r1.fq SRR2838702_R1.fastq.gz
            compress subsample-r2.fq SRR2838702_R2.fastq.gz
        fi

        if has_ont_reads; then
            compress subsample-ont.fq SRR2838702_ONT.fastq.gz
        fi

        if is_illumina_se; then
            compress subsample-r1.fq SRR2838702_SE.fastq.gz
        fi

        if [ "false" == "false" ]; then
            # Remove remaining intermediate FASTQ files
            rm *.fq
        fi
    fi
fi

#==========================================================================================
# Calculate quality stats before and after QC, optionally include plots (Default: true)
#==========================================================================================
if [ "${ERROR}" -eq "0" ]; then
    # fastq-scan
    if is_paired; then
        fastq_scan staging/r1/SRR2838702_R1.fastq.gz supplemental/SRR2838702_R1-original.json
        fastq_scan staging/r2/SRR2838702_R2.fastq.gz supplemental/SRR2838702_R2-original.json
        fastq_scan SRR2838702_R1.fastq.gz supplemental/SRR2838702_R1-final.json
        fastq_scan SRR2838702_R2.fastq.gz supplemental/SRR2838702_R2-final.json
    fi

    if has_ont_reads; then
        fastq_scan null supplemental/SRR2838702_ONT-original.json
        fastq_scan SRR2838702_ONT.fastq.gz supplemental/SRR2838702_ONT-final.json
    fi

    if is_illumina_se; then
        fastq_scan staging/r1/SRR2838702_R1.fastq.gz supplemental/SRR2838702_SE-original.json
        fastq_scan SRR2838702_SE.fastq.gz supplemental/SRR2838702_SE-final.json
    fi

    # FastQC and NanoPlot
    if [[ "false" == "false" ]]; then
        if has_ont_reads; then
            run_nanoplot "null" "original"
            run_nanoplot "SRR2838702_ONT.fastq.gz" "final"
        fi

        if has_illumina_reads; then
            # Fix issue on HPC when /tmp not writable
            mkdir -p ./tmp
            fastqs=()
            if is_paired; then
                ln -s staging/r1/SRR2838702_R1.fastq.gz SRR2838702_R1-original.fastq.gz
                ln -s staging/r2/SRR2838702_R2.fastq.gz SRR2838702_R2-original.fastq.gz
                ln -s SRR2838702_R1.fastq.gz SRR2838702_R1-final.fastq.gz
                ln -s SRR2838702_R2.fastq.gz SRR2838702_R2-final.fastq.gz
                fastqs=(SRR2838702_R1-original.fastq.gz SRR2838702_R2-original.fastq.gz SRR2838702_R1-final.fastq.gz SRR2838702_R2-final.fastq.gz)
            else
                ln -s staging/r1/SRR2838702_R1.fastq.gz SRR2838702_SE-original.fastq.gz
                ln -s SRR2838702_SE.fastq.gz SRR2838702_SE-final.fastq.gz
                fastqs=(SRR2838702_SE-original.fastq.gz SRR2838702_SE-final.fastq.gz)
            fi

            # Run FastQC
            fastqc                     --noextract                     --dir ./tmp                     -f fastq                     -t 4 "${fastqs[@]}"
            mv *_fastqc.html *_fastqc.zip supplemental/
            rm -rf tmp/ *-original.fastq.gz *-final.fastq.gz
        fi
    fi
fi

#==========================================================================================
# Final QC check
#
# This step ensures that the final reads meet minimum requirements. Users can skip this
# step with --skip_fastq_check
#==========================================================================================
if [ "false" == "false" ]; then
    # Only check for errors if we haven't already found them
    if [ "${ERROR}" -eq "0" ]; then
        # Check minimum coverage is met
        if [ "paired-end" == "hybrid" ]; then
            # Hybrid assembly, base checks on Illumina Reads
            gzip -cd SRR2838702_R1.fastq.gz SRR2838702_R2.fastq.gz | fastq-scan -g 358242 > temp.json
        elif [ "paired-end" == "short_polish" ]; then
            # Short read polishing, base checks on Nanopore Reads
            fastq_scan SRR2838702_ONT.fastq.gz temp.json
        else
            # base checks on which ever reads are available
            gzip -cd *.fastq.gz | fastq-scan -g 358242 > temp.json
        fi
        FINAL_BP=$(grep "total_bp" temp.json | sed -r 's/.*:[ ]*([0-9]+),/\1/')
        MIN_COVERAGE=$(( 10*358242 ))
        if [ ${FINAL_BP} -lt ${MIN_COVERAGE} ]; then
            echo "After QC, SRR2838702 FASTQ(s) contain ${FINAL_BP} total basepairs. This does
                    not exceed the required minimum ${MIN_COVERAGE} bp [10x coverage]. Further analysis
                    is discontinued." | \
            sed 's/^\s*//' > SRR2838702-low-sequence-coverage-error.txt
            ERROR=2
        fi
        rm temp.json

        # Final QC checks
        OPTS="--sample SRR2838702 --min_basepairs 2241820 --min_reads 7472 --min_proportion 0.5 --runtype paired-end"
        if is_paired; then
            # Check paired-end reads have same read counts
            if ! reformat.sh in1=SRR2838702_R1.fastq.gz in2=SRR2838702_R2.fastq.gz qin=auto out=/dev/null 2> SRR2838702-paired-end-error.txt; then
                echo "SRR2838702 FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \
                sed 's/^\s*//' >> SRR2838702-paired-end-error.txt
                ERROR=2
            else
                rm -f SRR2838702-paired-end-error.txt
            fi
            bactopia-check-fastqs --fq1 supplemental/SRR2838702_R1-final.json  --fq2 supplemental/SRR2838702_R2-final.json ${OPTS}
        fi

        if has_ont_reads; then
            bactopia-check-fastqs --fq1 supplemental/SRR2838702_ONT-final.json ${OPTS}
        fi

        if is_illumina_se; then
            bactopia-check-fastqs --fq1 supplemental/SRR2838702_SE-final.json ${OPTS}
        fi

        if ls *-error.txt 1>/dev/null 2>&1; then
            ERROR=2
        fi

    fi
fi

if is_assembly; then
    touch supplemental/reads-simulated-from-assembly.txt
fi

#==========================================================================================
# Handle Errors
#
# ERROR=0 : QC process successful
# ERROR=1 : Original reads could not be processed (e.g. empty after trimming)
# ERROR=2 : Final QC checks failed (e.g. insufficient reads after QC)
#==========================================================================================
if [ "${ERROR}" -eq "1" ]; then
    # Use original reads since we couldn't process them
    if is_paired; then
        cp staging/r1/SRR2838702_R1.fastq.gz supplemental/SRR2838702_R1.error-fastq.gz
        cp staging/r2/SRR2838702_R2.fastq.gz supplemental/SRR2838702_R2.error-fastq.gz

        if [ -s repair-singles.fq ]; then
            compress repair-singles.fq supplemental/SRR2838702.singles.fastq.gz
        fi
    fi

    if has_ont_reads; then
        cp null supplemental/SRR2838702_ONT.error-fastq.gz
    fi

    if is_illumina_se; then
        cp staging/r1/SRR2838702_R1.fastq.gz supplemental/SRR2838702_SE.error-fastq.gz
    fi
elif [ "${ERROR}" -eq "2" ]; then
    # Use the final processed reads which failed QC checks
    if is_paired; then
        mv SRR2838702_R1.fastq.gz supplemental/SRR2838702_R1.error-fastq.gz
        mv SRR2838702_R2.fastq.gz supplemental/SRR2838702_R2.error-fastq.gz

        if [ -s repair-singles.fq ]; then
            compress repair-singles.fq supplemental/SRR2838702.singles-error.fastq.gz
        fi
    fi

    if has_ont_reads; then
        mv SRR2838702_ONT.fastq.gz supplemental/SRR2838702_ONT.error-fastq.gz
    fi

    if is_illumina_se; then
        mv SRR2838702_SE.fastq.gz supplemental/SRR2838702_SE.error-fastq.gz
    fi
fi

# Cleanup

cat <<-END_VERSIONS > versions.yml
"QC:QC_MODULE":
    bbduk: $(echo $(bbduk.sh --version 2>&1) | sed 's/^.*BBTools version //;s/ .*$//')
    fastp: $(echo $(fastp --version 2>&1) | sed -e "s/fastp //g")
    fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
    lighter: $(echo $(lighter -v 2>&1) | sed 's/Lighter v//')
    nanoq: $(echo $(nanoq --version 2>&1) | sed 's/nanoq //')
    pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    porechop: $(echo $(porechop --version 2>&1))
    rasusa: $(echo $(rasusa --version 2>&1) | sed 's/rasusa //')
END_VERSIONS

if [[ "false" == "false" ]]; then
cat <<-END_EXTRA >> versions.yml
"QC:QC_MODULE":
    fastqc: $(echo $(fastqc --version 2>&1) | sed 's/^.*FastQC v//')
    nanoplot: $(echo $(NanoPlot -v 2>&1) | sed 's/.*NanoPlot //')
END_EXTRA
fi
