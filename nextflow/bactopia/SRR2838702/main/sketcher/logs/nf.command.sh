#!/bin/bash -ue
if [ "true" == "true" ]; then
    xz -c -d mash-refseq88.k21.msh.xz > ./mash-refseq88.k21.msh
fi

gzip -cd SRR2838702.fna.gz | mash sketch -o SRR2838702-k21 -k 21 -s 10000 -I SRR2838702 -
gzip -cd SRR2838702.fna.gz | mash sketch -o SRR2838702-k31 -k 31 -s 10000 -I SRR2838702 -
sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=10000 --merge SRR2838702 -o SRR2838702.sig SRR2838702.fna.gz

# Mash Screen
echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/	/g' > SRR2838702-mash-refseq88-k21.txt
gzip -cd SRR2838702.fna.gz | mash screen -w -i 0.8 -p 4 mash-refseq88.k21.msh - | sort -gr >> SRR2838702-mash-refseq88-k21.txt

# Sourmash classify
sourmash lca classify --query SRR2838702.sig --db gtdb-rs207.genomic-reps.dna.k31.lca.json.gz > SRR2838702-sourmash-gtdb-rs207-k31.txt

# Cleanup
if [ "true" == "true" ]; then
    rm -rf mash-refseq88.k21.msh
fi

# Capture versions
cat <<-END_VERSIONS > versions.yml
"SKETCHER:SKETCHER_MODULE":
    mash: $(echo $(mash 2>&1) | sed 's/^.*Mash version //;s/ .*$//')
    sourmash: $(echo $(sourmash --version 2>&1) | sed 's/.*sourmash //;')
END_VERSIONS
