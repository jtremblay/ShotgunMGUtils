# ShotgunMGUtils
Utilities to organize and query results from the ShotgunMG pipeline.

As sequencing output becomes increasingly larger, so are end results coming out of a bioinformatic pipeline. This utility was written to convert key result files coming out of the ShotgunMG pipeline. Briefly, it takes the annotation matrix, gene abundance matrix and metadata to create one single hdf5 file that can then be repeatedly (and rapidly) queried for downstream analysis. The motivation for this project was that it is not practical to load a gene abundance matrix of 15 GB into R as it consumes lots of RAM. In a hdf5 file format, however, RAM consumption is not an issue.

First, we have to build the hdf5 database:
```
~/build/ShotgunMGUtils$ ./shotgunmg_utils.py build \
    -m ./data/mapping_file.tsv \
    -g ./data/gene_abundance.tsv \
    -a ./data/annotations.tsv \
    -d ./db/contigs.h5 \
    -j ./db/index_ann.pckl \
    -k ./db/index_gene_abun.pckl
```
You should see something like this in output:
```
build
[DEBUG] Number of genes: 63592
[DEBUG] Number of samples: 6
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 63592/63592 [00:00<00:00, 89296.62it/s]
['Microbial_Mock_Community_Staggered_A', 'Microbial_Mock_Community_Even_C', 'Microbial_Mock_Community_Even_A', 'Microbial_Mock_Community_Even_B', 'Microbial_Mock_Community_Staggered_C', 'Microbial_Mock_Community_Staggered_B']
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 63592/63592 [00:00<00:00, 116242.68it/s]
['Treatment', 'Rep']
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 15401.36it/s]
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 63592/63592 [00:00<00:00, 166528.17it/s]
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 63592/63592 [00:00<00:00, 351814.35it/s]
Closing remaining open files:./db/contigs.h5...done
```
Then, we want to query the hdf5 database with something like this:
```
./shotgunmg_utils.py query \
    -d ./db/contigs.h5 \
    -j ./db/index_ann.pckl \
    -k ./db/index_gene_abun.pckl \
    --variable Treatment --variable-values Staggered:Even \
    --taxon genus --taxon-value Pseudomonas \
    --show_KO
```
On the standard error output, we get the following:
```
Showing aggregated gene counts 
 for (Treatment) Staggered and Even
 for (genus) Pseudomonas
 A total of 76 genes where found matching these criteria.
     Even |     19345.67 █████████████████████████
Staggered |      9580.00 ████████████▍

Occurence of most 20 abundant KEGG orthologs
b'K03406' |         6.00 █████████████████████████
b'K15125' |         4.00 ████████████████▋
b'K11904' |         2.00 ████████▎
b'K11903' |         2.00 ████████▎
b'K00138' |         1.00 ████▏
b'K03547' |         1.00 ████▏
b'K03546' |         1.00 ████▏
b'K03581' |         1.00 ████▏
b'K03582' |         1.00 ████▏
b'K03583' |         1.00 ████▏
b'K07240' |         1.00 ████▏
b'K16322' |         1.00 ████▏
b'K02278' |         1.00 ████▏
b'K12511' |         1.00 ████▏
b'K12510' |         1.00 ████▏
b'K02283' |         1.00 ████▏
b'K02282' |         1.00 ████▏
b'K02280' |         1.00 ████▏
b'K02279' |         1.00 ████▏
b'K02651' |         1.00 ████▏
Closing remaining open files:./db/contigs.h5...done
```

The resulting ```./db/contigs.h5``` file can be used in R with the appropriate libraries to generate more sophisticated analyses without being limited in RAM.



