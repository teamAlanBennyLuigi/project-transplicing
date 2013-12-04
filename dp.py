import os
import sys
import time

# Check for valid number of arguments
if (len(sys.argv) == 3):
    seqfile = sys.argv[1]
    reffile = sys.argv[2]
else:
    print 'Default sequence file and reference sequence used'
    print 'Usage: python findtrassplicingdp.py <seqfile.fasta> <reffile.fasta>'
    seqfile = "test.fasta"
    reffile = "sequence6992.fasta"

# Create output filename with unique id
uniqueid_t = time.strftime("%y-%m-%d_%H%M")
outputfile = "match_dp_%s" % (uniqueid_t)

# Run module and save output to text
os.system("python -u getrawtranscriptassignment.py %s %s >> %s.txt"
          % (seqfile, reffile, outputfile) )
print '*** Completed. ***\n File saved to:\n %s.txt' % (outputfile)
