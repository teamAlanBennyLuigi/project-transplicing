from Bio import SeqIO
import sys

# Check for valid number of arguments
if (len(sys.argv) == 3):
    infastq = sys.argv[1]
    outfasta = sys.argv[2]
else:
    print 'Default input sequence used'
    print 'Usage: python fastqtofasta.py <input_seq.fastq> <output_seq.fasta>'
    infastq = "IonXpress_001_R_2012_08_24_14_21_57_user_SN1-11-LopezRNAseq082412_SN1-11-LopezRNAseq082412LowStr.fastq"
    outfasta = "test.fasta"

# Concert .fastq to .fasta
count = SeqIO.convert(infastq, "fastq", outfasta, "fasta")
print "Converted %i records" % count
