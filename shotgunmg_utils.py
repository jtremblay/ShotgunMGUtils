#!/usr/bin/env python
 
"""Software to build and access shotgun metagenomic data (i.e. results of ShotgunMG pipeline) in a hdf5 database.

Julien Tremblay - jtremblay514@gmail.com
"""
import os
import sys
import argparse
import re
import h5py
import numpy as np 

cwd = os.path.dirname(__file__)
sys.path.insert(0, cwd + '')
#print(sys.path)
from lib.shotgunmg_db import *
  
def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    
    parser_build = subparsers.add_parser('build')
    parser_build.add_argument('-t', '--type', help="Type of database to build (default: shotgunmg)", choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_build.add_argument('-d', '--database-file', help="Database file", type=argparse.FileType('w'))
    parser_build.add_argument('-a', '--annotations-file', help="Annotations file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-g', '--gene-abundance-matrix-file', help="Gene abundance matrix file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-m', '--mapping-file', help="Mapping file used as input in the ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-j', '--index-annotations-file', help="Link pickle file for linking gene annotation to row", type=argparse.FileType('w'))
    parser_build.add_argument('-k', '--index-gene-abundance-file', help="Link pickle file for linking gene if to row", type=argparse.FileType('w'))
    
    parser_query = subparsers.add_parser('query')
    parser_query.add_argument('-t', '--type', help='Type of database to build (default: shotgunmg)', choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_query.add_argument('-d', '--database-file', help="Database file", type=argparse.FileType('r'))
    parser_query.add_argument('-o', '--outfile', help="Output file, for <query>", type=argparse.FileType('r'))
    parser_query.add_argument('-j', '--index-annotations-file', help="Link pickle file for linking gene annotation to row", type=argparse.FileType('r'))
    parser_query.add_argument('-k', '--index-gene-abundance-file', help="Link pickle file for linking gene if to row", type=argparse.FileType('r'))
    parser_query.add_argument('-v', '--variable', help="Variable (i.e. columns) to select in the mapping file", type=str)
    parser_query.add_argument('-z', '--variable-values', help="Two values, separated by a ':' character, that are present in the selected variable (arg --variable) from the mapping file to select", type=str)
    parser_query.add_argument('-l', '--taxon', help="Taxon level to select from the annotations", choices=["kingdom", "phylum", "class", "order", "family", "genus", "species"])
    parser_query.add_argument('-x', '--taxon-value', help="Taxon value to select from the taxon level (arg -l) (e.g. Lactobacillus)", type=str)
    parser_query.add_argument('-s', '--show_KO', help="Show histogram of most abundant KO satisfying -v, -z, -l, -x conditions.", action='store_true')
    
    parser_build = subparsers.add_parser('convert_abundance_to_json', help="Will convert the output of the gene_abundance file to json format and dump it in standard output")
    parser_build.add_argument('-t', '--type', help="Type of database to build (default: shotgunmg)", choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_build.add_argument('-d', '--database-file', required=False, help="Database file", type=argparse.FileType('w'))
    parser_build.add_argument('-a', '--annotations-file', required=False, help="Annotations file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-g', '--gene-abundance-matrix-file', required=True, help="Gene abundance matrix file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-m', '--mapping-file', required=False, help="Mapping file used as input in the ShotgunMG pipeline output", type=argparse.FileType('r'))
    
    parser_build = subparsers.add_parser('convert_annotations_to_json', help="Will convert the output of the annotations file to json format and dump it in standard output")
    parser_build.add_argument('-t', '--type', help="Type of database to build (default: shotgunmg)", choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_build.add_argument('-d', '--database-file', required=False, help="Database file", type=argparse.FileType('w'))
    parser_build.add_argument('-a', '--annotations-file', required=True, help="Annotations file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-g', '--gene-abundance-matrix-file', required=False, help="Gene abundance matrix file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-m', '--mapping-file', required=False, help="Mapping file used as input in the ShotgunMG pipeline output", type=argparse.FileType('r'))
    
    args = parser.parse_args(arguments)
    #exit(0);   
    if args.database_file is None and args.command != "convert_abundance_to_json" and args.command != "convert_annotations_to_json":
        raise ValueError('--database-file needed')
    
    if args.command == 'build':
        print("[INFO] build", file=sys.stderr)
        shotgun_mg = ShotgunMG(args.annotations_file.name, 
                               args.gene_abundance_matrix_file.name, 
                               args.mapping_file.name,
                               args.database_file.name,
                               write = 1,
                               dsetname = "Mock_community") 
        
        shotgun_mg.populate_annotations()
        shotgun_mg.populate_abundance()
        shotgun_mg.populate_mapping()
        shotgun_mg.build_indexes(args.index_annotations_file.name,
                                 args.index_gene_abundance_file.name)

    elif args.command == 'query':
        print("[INFO] query", file=sys.stderr)
        shotgun_mg = ShotgunMG(h5db_file = args.database_file.name,
                               write = 0)

        """Load indexes"""
        shotgun_mg.load_indexes(args.index_annotations_file.name,
                                 args.index_gene_abundance_file.name)

        """Store -v, -z, -l -x args in variales"""
        variable = args.variable
        variable_values = args.variable_values
        variable_value1 = args.variable_values.split(":")[0]
        variable_value2 = args.variable_values.split(":")[1]
        taxon = args.taxon
        taxon_value = args.taxon_value
        # For instance. Pull our all genes that are 2x fold change between
        # Staggered vs Even for Pseudomonas genus. Display results as average of both conditions
        # or display counts per column. maybe do ascii plots as well.
        shotgun_mg.query_v3(None, [variable, variable_value1, variable_value2], [taxon, taxon_value], [None, None], breakdown_by_KO=args.show_KO)
        shotgun_mg.barchart()
    
    if args.command == 'convert_abundance_to_json':
        print("[INFO] convert_abundance_to_json", file=sys.stderr)
        shotgun_mg = ShotgunMG("",#args.annotations_file.name, 
                               args.gene_abundance_matrix_file.name, 
                               "",#args.mapping_file.name,
                               "",
                               write = 2,
                               dsetname = "Mock_community") 
        
        shotgun_mg.generate_json_abundance()
    
    if args.command == 'convert_annotations_to_json':
        print("[INFO] convert_annotations_to_json", file=sys.stderr)
        shotgun_mg = ShotgunMG(args.annotations_file.name, 
                               "",#args.gene_abundance_matrix_file.name, 
                               "",#args.mapping_file.name,
                               "",
                               write = 2,
                               dsetname = "Mock_community") 
        
        shotgun_mg.generate_json_annotations()
     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


