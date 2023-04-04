# Build hdf5 database
mkdir -p db

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
