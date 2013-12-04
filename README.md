712 Project
===========
02-712: Computational Methods for Biological Modeling and Simulation  
Project: Interchromosomal and Intrachromosomal Trans-splicing Detection

**Team members:**
+ Alan Shteyman
+ Benny Jacob
+ Luigi Leung

## Usage
Concert .fastq to .fasta  
Defaults converts IonXpress_001_[...].fastq to test1.fasta  
.fastq are Drosophila transcripts from Professor Lopez  

        python fastqtofasta.py input_seq.fastq output_seq.fasta
        (or for defaults)
        python fastaqtofasta.py


Detect Trans-splicing events  
Using Dynamic Programming method: getrawtranscriptassignment.py  
Defaults using reference transcriptome sequence6992.fasta from modENCODE  
(from the [Berkeley Drosophila Genome Project]
(http://www.fruitfly.org/sequence/release5genomic.shtml))

        Defaults:
        InputSequence (seqfile) = test1.fasta
        ReferenceSequence (reffile) = sequence6992.fasta

        python getrawtranscriptassignment.py seqfile.fasta reffile.fasta
        (or for defaults)
        python getrawtranscriptassignment.py

Detect Trans-splicing events  
Using BLAST method: TBA
