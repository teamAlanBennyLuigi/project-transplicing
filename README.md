712 Project
===========
_02-712: Computational Methods for Biological Modeling and Simulation_  
**Project: Interchromosomal and Intrachromosomal Trans-splicing Detection**

**Team members:**
+ Alan Shteyman
+ Benny Jacob
+ Luigi Leung

Usage
-----
### Convert `.fastq` to `.fasta`
Defaults converts `IonXpress_001_[...].fastq` to `test.fasta`  
`IonXpress_001_[...].fastq` are Drosophila transcripts from [Professor Lopez]
(http://www.cmu.edu/bio/faculty/lopez.html)  

        $ python fastqtofasta.py input_seq.fastq output_seq.fasta
        (or for defaults)
        $ python fastaqtofasta.py


### Detect Trans-splicing events: DP
**Using Dynamic Programming method:** `getrawtranscriptassignment.py`  
Defaults using `test.fasta` that was generated from the previous step and  
reference transcriptome `sequence6992.fasta` from modENCODE
([Berkeley Drosophila Genome Project]
(http://www.fruitfly.org/sequence/release5genomic.shtml))

        Defaults:
        InputSequence (seqfile) = test.fasta
        ReferenceSequence (reffile) = sequence6992.fasta

        $ python dp.py seqfile.fasta reffile.fasta
        (or for defaults)
        $ python dp.py

### Detect Trans-splicing events: BLAST
**Using BLAST method:** TBA
