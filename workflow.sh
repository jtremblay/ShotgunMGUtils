# Build hdf5 database
mkdir -p db

# Validate files and convert to json
./shotgunmg_utils.py convert_to_json \
    -m ./data/mapping_file.tsv \
    -g ./data/gene_abundance.tsv \
    -a ./data/annotations.tsv > db/output.json

./shotgunmg_utils.py build \
    -m ./data/mapping_file.tsv \
    -g ./data/gene_abundance.tsv \
    -a ./data/annotations.tsv \
    -d ./db/contigs.h5 \
    -j ./db/index_ann.pckl \
    -k ./db/index_gene_abun.pckl

# Query the newly created database.
./shotgunmg_utils.py query \
    -d ./db/contigs.h5 -j ./db/index_ann.pckl -k ./db/index_gene_abun.pckl \
    --variable Treatment --variable-values Staggered:Even \
    --taxon genus --taxon-value Pseudomonas --show_KO

# If you wish to convert annotations data to json format:
./shotgunmg_utils.py convert_annotations_to_json \
    -a ./data/annotations.tsv  > ./db/annotations.json

# and convert gene abundance data to json format:
./shotgunmg_utils.py convert_abundance_to_json \
    -g ./data/gene_abundance.tsv > ./db/annotations.json
