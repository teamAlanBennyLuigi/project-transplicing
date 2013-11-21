#!/usr/bin/env python3

import os
import sys

# This is another way that runs on my computer. Try on your's
cmd = "blastn -out hehe.xml -outfmt 5 -query test1.fasta -db db/refdb -evalue 0.001"
print(cmd)
os.system(cmd)

#===============================================================================
# A more flexible version for later
#===============================================================================
#def blastn_cline(input_fasta, db_path, output_xml):
#    cmd = "blastn -out %s.xml -outfmt 5 -query %s.fasta -db %s -evalue 0.001" % input_fasta, output_xml, db_path
#    os.system(cmd)
#
#def main():
#    input_fasta = "test1"
#    db_path = "db/refdb"
#    output_xml = "keke"
#
#    blastn_cline("test1", "db/refdb", "keke")
#
#if __name__ == "__main__":
#    main()


#===============================================================================
# Benny's
#===============================================================================
#from Bio.Blast.Applications import NcbiblastnCommandline
#blastn_cline = NcbiblastnCommandline(query=r'test1.fasta', db=r"db/refdb",evalue=0.001,outfmt=5, out="hehe.xml")
#print blastn_cline
#blastn_cline()

