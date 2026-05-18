# DMScall
Single mutation variant calling tool for DMS sequencing data

#### What is this pipeline used for?
Deep mutational scanning (DMS) is a high-throughput screening method to comprehensively investigate the role of each amino acid residue in the protein. In this repository we showcase the analysis pipeline for the sequencing data acquired from DMS on STING (stimulator of interferon genes).

#### Step 1: Search and translate (mycall1.py)
**INPUT**: <br>
**fixed_seq**: Uniform nucleotide sequence around 7 nt right before the segment in which mutagenesis is introduced. <br>
**aa_length**: Peptide length of the segment in which mutagenesis is introduced. <br>
**input_fastq**: Merged paired-end sequencing data in .fastq format.<br>
<br>
**OUTPUT**:<br>
**output_fasta**: Nucleotide sequencing data translated into peptide sequences in .txt format. <br>
