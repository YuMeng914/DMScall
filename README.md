# DMScall
Single mutation variant calling tool for DMS sequencing data

### What is this pipeline used for?
Deep mutational scanning (DMS) is a high-throughput screening method to comprehensively investigate the role of each amino acid residue in the protein. In this repository we showcase the analysis pipeline for the sequencing data acquired from DMS on STING (stimulator of interferon genes).

### Step 1: Search and translation (mycall1.py)
**INPUT**: <br>
**fixed_seq**: Uniform nucleotide sequence around 7 nt right before the segment in which mutagenesis is introduced. <br>
**aa_length**: Peptide length of the segment in which mutagenesis is introduced. <br>
**input_fastq**: Merged paired-end sequencing data in .fastq format.<br>
<br>
**OUTPUT**:<br>
**output_fasta**: Nucleotide sequencing data translated into peptide sequences in .txt format. <br>

### Step 2: Alignment and mutation counting (mycall2.py)
**INPUT**:<br>
**input_file**: Translated sequencing data, acquired from Step 1. <br>
**reference_seq**: Wild-type peptide sequence of the segment in which mutagenesis is introduced. <br>
<br>
**OUTPUT**:<br>
**mutation_counts**: .txt file registering the amino-acid mutation counts across all translated sequencing reads relative to the wild-type reference sequence. <br>
**output_file**: Translated reads with only one amino-acid mutation compared to the wild-type sequence in .txt format. <br>

### Step 3: Variant calling (mycall3.py)
**INPUT**:<br>
**input_file**: Translated reads with only one mutation, acquired from Step 2. <br>
**reference_seq**: Same as in Step 2. <br>
<br>
**OUTPUT**:<br>
**mutation_matrix**: Variant calling matrix in .csv format. <br>

### Step 4: Enrichment analysis and visualization (mycall4_for_wt_normalization.py)
**INPUT**:<br>
**input_file_ctrl**: Variant calling matrix of ctrl group of cells, acquired in Step 3. <br>
**input_file_high**: Variant calling matrix of selected group of cells, acquired in Step 3. <br>
**ctrl_readcounts**: **mutation_counts** file of ctrl group of cells, acquired in Step 2. <br>
**high_readcounts**: **mutation_counts** file of selected group of cells, acquired in Step 2. <br>
**wt_sequence**: Wild-type peptide sequence of the segment in which mutagenesis is introduced, same as in Step 2 and Step 3. <br>
<br>
**PARAMETER**:<br>
**threshold**: Cut-off number of reads containing a specific single amino-acid mutation. <br>
<br>
**OUTPUT**:<br>
**savepath**: Heatmap visualizing enrichment scores of each mutation in .pdf format. <br>

