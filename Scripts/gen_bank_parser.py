# -*- coding: utf-8 -*-
"""
Created on Wed Jun 06 15:05:14 2018
Description: Parses the GenBank format file and extracts the CDS of the corresponding protein (hexon protein) and its protein id and writes it to an output file
@author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)
WARNING: This script assumes your input file ("RG_Test.gb") exists in the current working directory and it also writes the output to the current working directory itself.
"""

from Bio import SeqIO
import os

input_file = "RG_Test.gb"
output_file_name = "rg_test.txt"

if not os.path.exists(output_file_name):
    for rec in SeqIO.parse(input_file, "gb"):
        for feature in rec.features:
            for key, val in feature.qualifiers.items():
                if "hexon protein" in val:
                    with open(output_file_name, "a") as ofile:
                        ofile.write("Protein id: {0}\n{1}: {2}\n{3}\n\n".format(feature.qualifiers['protein_id'][0], key, val[0], feature.qualifiers['translation'][0]))
else:
    print ("The output file already seem to exist in the current workding directory {0}. Please change the name of the output file".format(os.getcwd()))                    
                    