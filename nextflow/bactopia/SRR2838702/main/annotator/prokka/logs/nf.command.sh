#!/bin/bash -ue
if [ "true" == "true" ]; then
    gzip -c -d staging/fna/SRR2838702.fna.gz > staging/fna/SRR2838702.fna
fi

if [ "false" == "true" ]; then
    export PROKKA_DBDIR=$(echo "$(which prokka | sed "s=/prokka==")/../db")
    env
    mkdir tmp_prokka/
    TMPDIR=tmp_prokka/ bactopia-prokka \
        --evalue 1e-09 --coverage 80 --centre Bactopia \
        --cpus 4 \
        --prefix SRR2838702 \
         \
        --locustag SRR2838702 \
        --proteins proteins.faa \
         \
        staging/fna/SRR2838702.fna
    rm -rf tmp_prokka/
else
    prokka \
        --evalue 1e-09 --coverage 80 --centre Bactopia \
        --cpus 4 \
        --prefix SRR2838702 \
         \
        --locustag SRR2838702 \
        --proteins proteins.faa \
         \
        staging/fna/SRR2838702.fna
fi

# Make blastdb of contigs, genes, proteins
mkdir blastdb
cat SRR2838702/SRR2838702.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for SRR2838702" -out blastdb/SRR2838702.fna
cat SRR2838702/SRR2838702.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for SRR2838702" -out blastdb/SRR2838702.ffn
cat SRR2838702/SRR2838702.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for SRR2838702" -out blastdb/SRR2838702.faa
tar -cvf - blastdb/ | gzip -c > SRR2838702/SRR2838702-blastdb.tar.gz

if [[ "false" == "false" ]]; then
    gzip SRR2838702/*.gff
    gzip SRR2838702/*.gbk
    gzip SRR2838702/*.fna
    gzip SRR2838702/*.faa
    gzip SRR2838702/*.ffn
    gzip SRR2838702/*.sqn
    gzip SRR2838702/*.fsa
    gzip SRR2838702/*.tbl
fi

mv SRR2838702/* ./

# Cleanup
if [ "true" == "true" ]; then
    rm -rf staging/fna/SRR2838702.fna
fi
rm -rf SRR2838702/ blastdb/ *.pdb *.pjs *.pot *.ptf *.pto

cat <<-END_VERSIONS > versions.yml
"PROKKA:PROKKA_MODULE":
    makeblastdb: $( echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    prokka: $( echo $(prokka --version 2>&1) | sed 's/^.*prokka //')
END_VERSIONS
