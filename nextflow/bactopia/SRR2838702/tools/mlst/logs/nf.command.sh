#!/bin/bash -ue
# Extract database
if [ "true" == "true" ]; then
    mkdir database/
    tar -xzf mlst.tar.gz -C database
    MLST_DB=$(find database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
else
    MLST_DB=$(find mlst.tar.gz/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
fi

mlst \
    --full \
    --label SRR2838702 \
    --threads 4 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --minid 95 --mincov 10 --minscore 50 \
    SRR2838702.fna.gz \
    > SRR2838702.tsv

if [[ -f "$MLST_DB/DB_VERSION" ]]; then
    DB_VERSION=$(cat $MLST_DB/DB_VERSION)
else
    DB_VERSION="custom database"
fi

# Cleanup
if [ "true" == "true" ]; then
    rm -rf database/
fi

cat <<-END_VERSIONS > versions.yml
"MLST:MLST_MODULE":
    mlst: $( echo $(mlst --version 2>&1) | sed 's/^.*mlst //' )
    database: $DB_VERSION
END_VERSIONS
