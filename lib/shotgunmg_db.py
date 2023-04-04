#!/usr/bin/env python

__author__     = "Julien Tremblay"
__copyright__  = "Copyright 2023" 
__credits__    = ["Julien Tremblay"]
__license__    = "GPL"
__version__    = "1.0.0"
__maintainer__ = "Julien Tremblay"
__email__      = "jtremblay514@gmail.com"
__status__     = "Release"


"""Contains code for the ShotgunMG uilities project.
"""
import sys
import os
import csv
import pickle
from tqdm import *
import h5py
import tables as tb
import numpy as np 
import statistics
import warnings
from typing import List, Dict, Optional
from pydantic import BaseModel, ValidationError, create_model
warnings.filterwarnings('ignore', category=tb.NaturalNameWarning)

def do_barchart(data, longest_label_length):
    max_value = data[max(data, key=data.get, default=None)]
    increment = max_value / 25
    
    for key in data:
        label = str(key)
        count = round(data[key], 2)
        # The ASCII block elements come in chunks of 8, so we work out how
        # many fractions of 8 we need.
        # https://en.wikipedia.org/wiki/Block_Elements
        bar_chunks, remainder = divmod(int(count * 8 / increment), 8)

        # First draw the full width chunks
        bar = '█' * bar_chunks

        # Then add the fractional part.  The Unicode code points for
        # block elements are (8/8), (7/8), (6/8), ... , so we need to
        # work backwards.
        if remainder > 0:
            bar += chr(ord('█') + (8 - remainder))

        # If the bar is empty, add a left one-eighth block
        bar = bar or  '▏'
        print(f'{label.rjust(longest_label_length)} | {count:#12.2f} {bar}')


def vprint(obj):
    sys.stderr.write("[DEBUG] " + str(obj) + "\n")


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


class ShotgunMG:

    def __init__(self, annotations_file="", abundance_file="", mapping_file="", h5db_file="", index_annotations_file="", index_gene_abundance_file="", write=0, number_of_genes=0, number_of_samples=0, dsetname=""):

        self._h5db_file = h5db_file
        self._annotations_file = annotations_file
        self._abundance_file = abundance_file
        self._dsetname = dsetname
        self._mapping_file = mapping_file
        self._index_annotations_file = index_annotations_file
        self._index_gene_abundance_file = index_gene_abundance_file
        self._index_annotations = None
        self._index_gene_abundance = None
        # Query
        self._query_mean_abundance_values_treatment_foreach_gene = None
        self._query_mean_abundance_values_control_foreach_gene = None
        self._query_variable = None
        self._query_treatment = None
        self._query_control = None
        self._query_tax_category = None
        self._query_tax_value = None
        self._query_function_category = None
        self._query_function_value = None
        self._fc = None
        self._query_gene_ids = None
        self._KOs = None
        self._COGs = None
        self._KOs_by_treatment = None
        self._KOs_by_control = None
        self._gene_id_to_KO = None
        self._KO_to_gene_id = None
        self._breakdown_by_KO = False

        if(write == 1):
            self._h5f_handle = tb.open_file(h5db_file, mode="w", title="")
        elif(write == 0):
            self._h5f_handle = tb.open_file(h5db_file, mode="r")
        elif(write == 2):
            self._h5f_handle = ""
        else:
            raise Exception("write argument has to be = 1 or = 0 or = 2")
       
        if(annotations_file != ""):
            with open(annotations_file) as f:
                first_line = f.readline()
                for line in f:
                    number_of_genes = number_of_genes + 1
                f.close()
            self._number_of_genes = number_of_genes
            sys.stderr.write("[DEBUG] Number of genes: " + str(self._number_of_genes) + "\n")

        elif(abundance_file != ""):
            with open(abundance_file) as f:
                first_line = f.readline()
                for line in f:
                    number_of_genes = number_of_genes + 1
                f.close()
            self._number_of_genes = number_of_genes
            sys.stderr.write("[DEBUG] Number of genes: " + str(self._number_of_genes) + "\n")

        if(mapping_file != ""):
            with open(mapping_file) as f:
                first_line = f.readline()
                for line in f:
                    number_of_samples = number_of_samples + 1
                f.close()
            self._number_of_samples = number_of_samples
            sys.stderr.write("[DEBUG] Number of samples: " + str(self._number_of_samples) + "\n")

    @property
    def number_of_genes(self):
        return self._number_of_genes
    
    @property
    def annotations_file(self):
        return self._annotations_file

    @annotations_file.setter
    def annotations_file(self, value):
        self._annotations_file = value
    
    @property
    def abundance_file(self):
        return self._abundance_file
    
    @property
    def get_mean_abundance_values_control_foreach_gene(self):
        return self._mean_query_abundance_values_control_foreach_gene
    
    @property
    def get_mean_abundance_values_treatment_foreach_gene(self):
        return self._mean_query_abundance_values_treatment_foreach_gene

    def populate_annotations(self):
        """Parses annotations.tsv output file from ShotgunMG pipeline and stores it into 
           a hdf5 database.

           Returns nothing at the moment.
        """

        annotations_file = self._annotations_file
        h5w = self._h5f_handle

        dataset_name = "annotations"
        
        class DataDescr(tb.IsDescription):
            #gene_id   = tb.StringCol(32, pos=0)
            gene_id   = tb.Int32Col(pos=0)
            contig_id = tb.StringCol(32, pos=1)
            KO        = tb.StringCol(32, pos=2)
            COG       = tb.StringCol(32, pos=3)
            PFAM      = tb.StringCol(32, pos=4)
            kingdom   = tb.StringCol(32, pos=5)
            phylum    = tb.StringCol(32, pos=6)
            Class     = tb.StringCol(32, pos=7)
            order     = tb.StringCol(32, pos=8)
            family    = tb.StringCol(32, pos=9)
            genus     = tb.StringCol(32, pos=10)
            species   = tb.StringCol(32, pos=11)

        with open(annotations_file) as f:
            # Skip header line
            first_line = f.readline()
            i = 1
            table = h5w.create_table("/", 'annotations', DataDescr, "Annotation table unindexed", expectedrows=self._number_of_genes, filters=tb.Filters(complevel=9, complib='blosc'))
            my_pointer = table.row

            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    
                    my_pointer['gene_id']   = int(str(line_list[1]).replace("gene_id_",""))
                    my_pointer['contig_id'] = line_list[0]
                    my_pointer['KO']        = line_list[3]
                    my_pointer['COG']       = line_list[9]
                    my_pointer['PFAM']      = line_list[16]
                    my_pointer['kingdom']   = line_list[26]
                    my_pointer['phylum']    = line_list[27]
                    my_pointer['Class']     = line_list[28]
                    my_pointer['order']     = line_list[29]
                    my_pointer['family']    = line_list[30]
                    my_pointer['genus']     = line_list[31]
                    my_pointer['species']   = line_list[32]
                        
                    my_pointer.append()
                    pbar.update(1)
                table.flush() 
            f.close()


    def populate_abundance(self):
        """Essentially stores a gene abundance matrix in the hdf5 file. PyTables version 
           Returns nothing at the moment.
        """
        
        abundance_file = self._abundance_file
        h5w = self._h5f_handle

        dataset_name = "gene_abundance"
        number_of_genes = self._number_of_genes
    
        # Create IsDescription structure only for gene id
        class DataDescr(tb.IsDescription):
            gene_id = tb.Int32Col(pos=0)

        with open(abundance_file) as f:
            # Skip header line
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line_list = first_line.split('\t')
            first_line_list.pop(0)
            print(first_line_list)
            i = 1
            # Here we dynamically populate the abundance tables.
            for el in first_line_list:
                DataDescr.columns[el] = tb.Int32Col(pos=i)
            
            table = h5w.create_table("/", 'gene_abundance', DataDescr, "Gene abundance matrix - unindexed", expectedrows=self._number_of_genes, filters=tb.Filters(complevel=9, complib='blosc'))
            my_pointer = table.row

            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    row_id = line_list.pop(0)
                    row_id = int(str(row_id).replace("gene_id_",""))
                    my_pointer['gene_id'] = row_id
                        
                    for j in range(len(line_list)):
                        my_pointer[first_line_list[j]] = int(line_list[j])
                    my_pointer.append()
                    pbar.update(1)
                table.flush()


    def populate_mapping(self):
        """Parses mapping_file that was used in input for the ShotgunMG pipeline and stores it into 
           a hdf5 database. ***Sample ID has to be the first column.***

           Returns nothing at the moment.
        """
        
        mapping_file = self._mapping_file
        h5w = self._h5f_handle
        dataset_name = "mapping"
        
        class DataDescr(tb.IsDescription):
            sample_id = tb.StringCol(256, pos=0)

        with open(mapping_file) as f:
            # Skip header line
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line_list = first_line.split('\t')
            first_line_list.pop(0)
            print(first_line_list)
            i = 1
            # Here we dynamically populate the abundance tables.
            for el in first_line_list:
                DataDescr.columns[el] = tb.StringCol(64, pos=i)
            
            table = h5w.create_table("/", 'mapping', DataDescr, "Metadata")
            my_pointer = table.row

            with tqdm(total=self._number_of_samples) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    row_id = line_list.pop(0)
                    my_pointer['sample_id'] = row_id
                        
                    for j in range(len(line_list)):
                        my_pointer[first_line_list[j]] = str(line_list[j])
                    my_pointer.append()
                    pbar.update(1)
                table.flush()


    def query_gene_id(self, gene_id):
        """Query a gene_id and get its annotations and abundance values.
        """
        h5r = self._h5f_handle
        rows = h5r.root.gene_abundance.read_where('(gene_id == ' + str(gene_id) + ') | (gene_id == 14701)')
        print(rows)


    def query_genus(self, genus):
        """Query a genus and get its annotations and abundance values.
        """

        h5r = self._h5f_handle
        rows = h5r.root.annotations.read_where('genus == b"' + genus + '"')
        print(rows)


    def query(self, fc, variables, taxonomy, function, breakdown_by_KO=False):
        """General query.
        """

        variable = variables[0]
        treatment = variables[1]
        control = variables[2]
        tax_category = taxonomy[0]
        tax_value = taxonomy[1]
        function_category = function[0]
        function_value = function[1]
        
        self._query_variable = variable
        self._query_treatment = treatment
        self._query_control = control
        self._query_tax_category = tax_category
        self._query_tax_value = tax_value
        self._query_function_category = function_category
        self._query_function_value = function_value
        self._fc = fc

        h5r = self._h5f_handle

        """Select gene_ids based on taxonomy or function
        """
        my_annotations_table = h5r.root.annotations
        condition = '(' + tax_category + ' == b"' + tax_value + '")'
        
        complete_gene_ids = []
        KOs = {}
        COGs = {}
        KO_to_gene_id = {}
        for record in my_annotations_table.where(condition):
            complete_gene_ids.append(record['gene_id'])
            if record['KO'] in KOs:
                KOs[record['KO']] += 1
                #KO_to_gene_id[record['KO']] += ";"str(record['gene_id']
            else:
                KOs[record['KO']] = 1
                KO_to_gene_id[record['KO']] = record['gene_id']
            
            if record['COG'] in COGs:
                COGs[record['COG']] += 1
            else:
                COGs[record['COG']] = 1
        
        self._KOs = KOs
        self._COGs = COGs
        self._query_gene_ids = complete_gene_ids

        """Narrow selection to selected samples
        """
        my_mapping_table = h5r.root.mapping
        
        treatment_variables = '(' + variable + ' == b"' + treatment + '")'
        treatment_samples = []
        for record in my_mapping_table.where(treatment_variables):
            treatment_samples.append(record['sample_id'])
        treatment_samples = [str(ts).replace("b'","").replace("'","") for ts in treatment_samples]
        
        control_variables = '(' + variable + ' == b"' + control + '")'
        control_samples = []
        for record in my_mapping_table.where(control_variables):
            control_samples.append(record['sample_id'])
        control_samples = [str(ts).replace("b'","").replace("'","") for ts in control_samples]

        """Once selected genes (rows) and selected samples (cols) is done, 
           construct query string and get values from the gene abundance table
           for both treament and control.
        """
        my_abundance_table = h5r.root.gene_abundance
        # Example: row_condition = '(gene_id == b"gene_id_31516")'
        # Query by batches of 100, because otherwise request string will be too long.
        n = 31
        x = list(divide_chunks(complete_gene_ids, n))
        my_treatment_means = []
        my_treatment_stdevs = []
        my_control_means = []
        my_control_stdevs = []
        for gene_ids in x:
            row_condition = "(gene_id == "
            row_condition += " ) | (gene_id == ".join([str(gene_id) for gene_id in gene_ids])
            row_condition += ")"
           
            for row in my_abundance_table.where(row_condition, condvars=None):
                replicate_values = []
                for sample_name in treatment_samples:
                    #sys.stderr.write(sample_name + " " + str(row[sample_name]) + "\n")
                    replicate_values.append(row[sample_name])
                my_treatment_means.append(round(statistics.mean(replicate_values), 2))
                my_treatment_stdevs.append(round(statistics.stdev(replicate_values),2 ))
                
                replicate_values = []
                for sample_name in control_samples:
                    #sys.stderr.write(sample_name + " " + str(row[sample_name]) + "\n")
                    replicate_values.append(row[sample_name])
                my_control_means.append(round(statistics.mean(replicate_values), 2))
                my_control_stdevs.append(round(statistics.stdev(replicate_values),2 ))

        self._query_mean_abundance_values_treatment_foreach_gene = my_treatment_means
        self._query_mean_abundance_values_control_foreach_gene = my_control_means


    def query_v2(self, fc, variables, taxonomy, function, breakdown_by_KO=False):
        """General query.
        """

        variable = variables[0]
        treatment = variables[1]
        control = variables[2]
        tax_category = taxonomy[0]
        tax_value = taxonomy[1]
        function_category = function[0]
        function_value = function[1]
        
        self._query_variable = variable
        self._query_treatment = treatment
        self._query_control = control
        self._query_tax_category = tax_category
        self._query_tax_value = tax_value
        self._query_function_category = function_category
        self._query_function_value = function_value
        self._fc = fc

        h5r = self._h5f_handle

        """Select gene_ids based on taxonomy or function
        """
        my_annotations_table = h5r.root.annotations
        condition = '(' + tax_category + ' == b"' + tax_value + '")'
        vprint(condition)
        
        complete_gene_ids = []
        KOs = {}
        COGs = {}
        for record in my_annotations_table.where(condition):
            complete_gene_ids.append(record['gene_id'])
            if record['KO'] in KOs:
                KOs[record['KO']] += 1
            else:
                KOs[record['KO']] = 1
            
            if record['COG'] in COGs:
                COGs[record['COG']] += 1
            else:
                COGs[record['COG']] = 1
        
        self._KOs = KOs
        self._COGs = COGs
        self._query_gene_ids = complete_gene_ids

        """Narrow selection to selected samples
        """
        my_mapping_table = h5r.root.mapping
        
        treatment_variables = '(' + variable + ' == b"' + treatment + '")'
        treatment_samples = []
        for record in my_mapping_table.where(treatment_variables):
            treatment_samples.append(record['sample_id'])
        treatment_samples = [str(ts).replace("b'","").replace("'","") for ts in treatment_samples]
        
        control_variables = '(' + variable + ' == b"' + control + '")'
        control_samples = []
        for record in my_mapping_table.where(control_variables):
            control_samples.append(record['sample_id'])
        control_samples = [str(ts).replace("b'","").replace("'","") for ts in control_samples]

        """Once selected genes (rows) and selected samples (cols) is done, 
           construct query string and get values from the gene abundance table
           for both treament and control.
        """
        my_abundance_table = h5r.root.gene_abundance
        # Example: row_condition = '(gene_id == b"gene_id_31516")'
        # Query by batches of 100, because otherwise request string will be too long.
        n = 31
        x = list(divide_chunks(complete_gene_ids, n))
        my_treatment_means = []
        my_treatment_stdevs = []
        my_control_means = []
        my_control_stdevs = []
        my_gene_row_list = []
        
        
        for gene_ids in x:
            row_condition = "(gene_id == "
            row_condition += " ) | (gene_id == ".join([str(gene_id) for gene_id in gene_ids])
            row_condition += ")"

            gene_rows = my_abundance_table.get_where_list(row_condition, condvars=None)
            my_abundance_table_select = my_abundance_table.read_coordinates(gene_rows)
            
            replicate_values = []
            for sample_name in treatment_samples:
                replicate_values.append(my_abundance_table_select[sample_name])
            my_treatment_means.extend(np.mean(replicate_values, axis=0))
            my_treatment_stdevs.extend(np.std(replicate_values, axis=0))

            replicate_values = []
            for sample_name in control_samples:
                replicate_values.append(my_abundance_table_select[sample_name])
            
            my_control_means.extend(np.mean(replicate_values, axis=0))
            my_control_stdevs.extend(np.std(replicate_values, axis=0))
        
        self._query_mean_abundance_values_treatment_foreach_gene = my_treatment_means
        self._query_mean_abundance_values_control_foreach_gene = my_control_means


    def query_v3(self, fc, variables, taxonomy, function, breakdown_by_KO=False):
        """General query.
        """

        self._breakdown_by_KO = breakdown_by_KO
        variable = variables[0]
        treatment = variables[1]
        control = variables[2]
        tax_category = taxonomy[0]
        tax_value = taxonomy[1]
        function_category = function[0]
        function_value = function[1]
        
        self._query_variable = variable
        self._query_treatment = treatment
        self._query_control = control
        self._query_tax_category = tax_category
        self._query_tax_value = tax_value
        self._query_function_category = function_category
        self._query_function_value = function_value
        self._fc = fc
        
        index_annotations_file = self._index_annotations_file
        index_gene_abundance_file = self._index_gene_abundance_file
        index_annotations = self._index_annotations 
        index_gene_abundance = self._index_gene_abundance

        h5r = self._h5f_handle

        """Select gene_ids based on taxonomy or function
        """
        my_annotations_table = h5r.root.annotations
        condition = '(' + tax_category + ' == b"' + tax_value + '")'
        vprint(condition)
        
        complete_gene_ids = []
        KOs = {}
        COGs = {}
        for record in my_annotations_table.where(condition):
            complete_gene_ids.append(record['gene_id'])
            if record['KO'] in KOs:
                KOs[record['KO']] += 1
            else:
                KOs[record['KO']] = 1
            
            if record['COG'] in COGs:
                COGs[record['COG']] += 1
            else:
                COGs[record['COG']] = 1
        
        self._KOs = KOs
        self._COGs = COGs
        self._query_gene_ids = complete_gene_ids

        """Narrow selection to selected samples
        """
        my_mapping_table = h5r.root.mapping
        
        treatment_variables = '(' + variable + ' == b"' + treatment + '")'
        treatment_samples = []
        for record in my_mapping_table.where(treatment_variables):
            treatment_samples.append(record['sample_id'])
        treatment_samples = [str(ts).replace("b'","").replace("'","") for ts in treatment_samples]
        
        control_variables = '(' + variable + ' == b"' + control + '")'
        control_samples = []
        for record in my_mapping_table.where(control_variables):
            #print(record['sample_id'])
            control_samples.append(record['sample_id'])
        control_samples = [str(ts).replace("b'","").replace("'","") for ts in control_samples]

        """Once selected genes (rows) and selected samples (cols) is done, 
           construct query string and get values from the gene abundance table
           for both treament and control.
        """
        my_abundance_table = h5r.root.gene_abundance
        # Example: row_condition = '(gene_id == b"gene_id_31516")'
        # Query by batches of 100, because otherwise request string will be too long.
        n = 31
        x = list(divide_chunks(complete_gene_ids, n))
        my_treatment_means = []
        my_treatment_stdevs = []
        my_control_means = []
        my_control_stdevs = []
        my_gene_row_list = []
       
        gene_rows = []
        for gene_id in complete_gene_ids:
            gene_rows.append(self._index_gene_abundance[gene_id])
        
        my_abundance_table_select = my_abundance_table.read_coordinates(gene_rows)
        
        replicate_values = []
        for sample_name in treatment_samples:
            replicate_values.append(my_abundance_table_select[sample_name])
        
        my_treatment_means.extend(np.mean(replicate_values, axis=0))
        my_treatment_stdevs.extend(np.std(replicate_values, axis=0))

        replicate_values = []
        for sample_name in control_samples:
            replicate_values.append(my_abundance_table_select[sample_name])
        
        my_control_means.extend(np.mean(replicate_values, axis=0))
        my_control_stdevs.extend(np.std(replicate_values, axis=0))
        
        self._query_mean_abundance_values_treatment_foreach_gene = my_treatment_means
        self._query_mean_abundance_values_control_foreach_gene = my_control_means


    def barchart(self):
        
        data_control = np.sum(self._query_mean_abundance_values_control_foreach_gene)
        data_treatment = np.sum(self._query_mean_abundance_values_treatment_foreach_gene)
        data = {}
        data[self._query_control] = data_control
        data[self._query_treatment] = data_treatment
        #print(data)

        longest_label_length = max(len(self._query_treatment), len(self._query_control))
        print(longest_label_length)
        
        print("Showing aggregated gene counts ")
        print(" for (" + self._query_variable + ") " + self._query_treatment + " and " + self._query_control)
        print(" for (" + self._query_tax_category + ") " + self._query_tax_value)
        print(" A total of " + str(len(self._query_gene_ids)) + " genes where found matching these criteria.")
        do_barchart(data, longest_label_length)
        print("")

        if(self._breakdown_by_KO == True):
            KOs = self._KOs
            KOs.pop(b'NULL', None)
            KOs = dict(sorted(KOs.items(), key=lambda x: x[1], reverse=True))
            my_keys = list(KOs.keys())
            my_cutoff = 0
            if len(my_keys) >= 19 :
                my_cutoff = 20
            else:
                my_cutoff = len(my_keys)
            my_keys = my_keys[0:my_cutoff]
            KOs2 = {}
            for key in my_keys:
                KOs2[key] = KOs[key]
                
            longest_label_length = max(len(label) for label in KOs2)
            print("Occurence of most " + str(my_cutoff) + " abundant KEGG orthologs")
            do_barchart(KOs2, longest_label_length)


    def dump_mapping(self):
        """Dump mapping file.
        """
        h5r = self._h5f_handle
        table = h5r.root.mapping
        my_colnames = table.colnames
        my_dict = {}
        for i in range(len(my_colnames)):
            my_dict[str(i)] = my_colnames[i]

        print(my_dict)
        i = 0
        for r in table.iterrows():
            curr_header_name = my_dict[str(i)]
            print(r[curr_header_name], end='   ')
            
            i = i + 1
            if(i == len(my_colnames)):
               print("")
               i = 0

    def get_sample_ids_from_mapping(self, variable, value):
        """Dump mapping file.
        """
        h5r = self._h5f_handle
        table = h5r.root.mapping
        my_colnames = table.colnames
        my_dict = {}
        rows = h5r.root.mapping.read_where(variable + ' == b"' + value + '"')
        print("---")
        print(rows)


    """ Build index in pckl files
    """

    def build_indexes(self, index_annotations_file, index_gene_abundance_file):
        """indexes
           gene_id_n = row[x] in the h5 'annotations' and 'gene_abundance' datasets
        """
        index_annotations = dict()
        i = 0
        with tqdm(total=self._number_of_genes) as pbar:
            with open(self._annotations_file) as f:
                first_line = f.readline()
                for line in f:
                    line_list = line.split('\t')
                    gene_id = int(str(line_list[1]).replace("gene_id_",""))
                    index_annotations[gene_id] = i
                    pbar.update(1)
                    i = i + 1
            f.close()

        with open(index_annotations_file, 'wb') as f:
            pickle.dump(index_annotations, f)

        index_gene_abundance = dict()
        i = 0
        with tqdm(total=self._number_of_genes) as pbar:
            with open(self._abundance_file) as f:
                first_line = f.readline()
                for line in f:
                    line_list = line.split('\t')
                    gene_id = int(str(line_list[0]).replace("gene_id_",""))
                    index_gene_abundance[gene_id] = i
                    pbar.update(1)
                    i = i + 1
            f.close()

        with open(index_gene_abundance_file, 'wb') as f:
            pickle.dump(index_gene_abundance, f)

    def load_indexes(self, index_annotations_file, index_gene_abundance_file):
        """indexes
           gene_id_n = row[x] in the h5 'annotations' and 'gene_abundance' datasets
        """
        index_annotations = dict()
        with open(index_annotations_file, 'rb') as f:
            self._index_annotations = pickle.load(f)
        f.close()

        index_gene_abundance = dict()
        with open(index_gene_abundance_file, 'rb') as f: 
            self._index_gene_abundance = pickle.load(f)
        f.close()

    def generate_json_annotations(self):
        """Generating json output representing abundance and annotations files. Validate with pydantic.
           Will write json output to standard output in terminal. Returns nothing at the moment.
        """
        annotations_file = self._annotations_file
        dataset_name = "annotations"

        class Annotation(BaseModel):
            gene_id   : int
            contig_id : str
            KO        : str
            COG       : str
            PFAM      : str
            kingdom   : str
            phylum    : str
            Class     : str
            order     : str
            family    : str
            genus     : str
            species   : str
       
        with open(annotations_file) as f:
            # Skip header line
            first_line = f.readline()
            i = 1
            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')

                    data = {
                        'gene_id'   : int(str(line_list[1]).replace("gene_id_","")),
                        'contig_id' : line_list[0],
                        'KO'        : line_list[3],
                        'COG'       : line_list[9],
                        'PFAM'      : line_list[16],
                        'kingdom'   : line_list[26],
                        'phylum'    : line_list[27],
                        'Class'     : line_list[28],
                        'order'     : line_list[29],
                        'family'    : line_list[30],
                        'genus'     : line_list[31],
                        'species'   : line_list[32]
                    }
                    
                    try:
                        annotation_row = Annotation(**data)
                        print(annotation_row.dict(), file=sys.stdout)
                    except ValidationError as e:
                        print("Exception as str:", file=sys.stderr)
                        print(e)
                        print("Exception as json:", file=sys.stderr)
                        print(e.json())

                    pbar.update(1)
            f.close()
        

    def generate_json_abundance(self):
        """Generating json output representing abundance and annotations files. Validate with pydantic.
           Will write json output to standard output in terminal. Returns nothing at the moment.
        """
        # Dynamic model created on runtime
        Abundance = create_model(
           'Abundance',
           gene_id=(str, ...),
           abundance=(dict[str, float], ...)
        )
        
        with open(self._abundance_file) as f:
            # Skip header line
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line_list = first_line.split('\t')
            first_line_list.pop(0)
            i = 1
            # Here we dynamically populate the abundance tables.
            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    row_id = line_list.pop(0)
                    row_id = int(str(row_id).replace("gene_id_",""))
                    #my_pointer['gene_id'] = row_id
                        
                    data = {}
                    data["gene_id"] = row_id
                    data["abundance"] = {}
                    for j in range(len(line_list)):
                        #data[row_id][first_line_list[j]] = int(line_list[j])
                        data["abundance"][first_line_list[j]] = float(line_list[j])
    
                    try:
                        abundance_row = Abundance.parse_obj(data)
                        print(abundance_row.dict(), file=sys.stdout)
                    except ValidationError as e:
                        print("Exception as str:", file=sys.stderr)
                        print(e)
                        print("Exception as json:", file=sys.stderr)
                        print(e.json())
                    pbar.update(1)

