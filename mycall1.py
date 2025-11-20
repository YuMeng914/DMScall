from Bio import SeqIO
from Bio.Seq import Seq

fixed_seq = "gggatcc" #seg1
fixed_seq = "caacgtg" #seg2
fixed_seq = "gtttgcc" #seg3

aa_length = 155 #seg1
aa_length = 115 #seg2
aa_length = 109 #seg3

nt_length = aa_length * 3

input_fastq = "MSQ0314 Plasmid/group_3.fastq"
output_fasta = "MSQ0314 Plasmid/plas3_translated.txt"

with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fastq, "fastq"):
        seq = str(record.seq).lower()

        # Find fixed sequence
        pos = seq.find(fixed_seq)
        if pos == -1:
            # fixed sequence not found, skip read
            continue

        # Start after the fixed sequence
        start = pos + len(fixed_seq)
        if start + nt_length > len(seq):
            # not enough sequence after fixed_seq, skip
            continue

        # Extract subsequence for translation
        coding_seq = seq[start : start + nt_length]

        # Translate nucleotide sequence
        protein_seq = Seq(coding_seq).translate(to_stop=False)

        # Write protein to fasta
        out_f.write(f"{protein_seq}\n")