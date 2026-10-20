#!/bin/bash -ue
# Create a file of files for csvtk
ls staging/csv/* | awk '{ print $1 }' > fofn.txt

csvtk \
    concat \
      \
    --num-cpus 4 \
    --tabs  \
    --out-tabs \
    --out-file meta.tsv \
    --infile-list fofn.txt

# Cleanup
rm -rf fofn.txt

cat <<-END_VERSIONS > versions.yml
"GATHER:CSVTK_CONCAT":
    csvtk: $(echo $( csvtk version | sed -e "s/csvtk v//g" ))
END_VERSIONS
