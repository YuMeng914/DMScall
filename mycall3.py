import numpy as np
import pandas as pd

aa_list = list("ACDEFGHIKLMNPQRSTVWY")

reference_seq = "MPHSSLHPSIPCPRGHGAQKAALVLLSACLVTLWGLGEPPEHTLRYLVLHLASLQLGLLLNGVCSLAEELRHIHSRYRGSYWRTVRACLGCPLRRGALLLLSIYFYYSLPNAVGPPFTWMLALLGLSQALNILLGLKGLAPAEISAVCEKGNFNV"#seg1
#reference_seq = "AHGLAWSYYIGYLRLILPELQARIRTYNQHYNNLLRGAVSQRLYILLPLDCGVPDNLSMADPNIRFLDKLPQQTGDRAGIKDRVYSNSIYELLENGQRAGTCVLEYATPLQTLFA" #seg2
#reference_seq = "MSQYSQAGFSREDRLEQAKLFCRTLEDILADAPESQNNCRLIAYQEPADDSSFSLSQEVLRHLRQEEKEEVTVGSLKTSAVPSTSTMSQEPELLISGMEKPLPLRTDFS" #seg3

mutation_matrix = pd.DataFrame(
    np.nan,
    index=range(len(reference_seq)),  # 0 to 114
    columns=aa_list
)

flag = 0

input_file = "MSQ0314 Plasmid/1_one_mut.txt"
with open(input_file, "r") as f:
    for line_num, line in enumerate(f, 1):
        seq = line.strip()
        if len(seq) != len(reference_seq):
            print(f"Warning: line {line_num} has length {len(seq)}")
            continue

        # Find mutation position
        for i, (ref_aa, mut_aa) in enumerate(zip(reference_seq, seq)):
            if ref_aa != mut_aa:
                if mut_aa not in aa_list:
                    #print(f"Warning: non-standard AA '{mut_aa}' at line {line_num}, position {i}")
                    flag += 1
                    continue
                # Initialize cell if it's NaN
                if pd.isna(mutation_matrix.at[i, mut_aa]):
                    mutation_matrix.at[i, mut_aa] = 1
                else:
                    mutation_matrix.at[i, mut_aa] += 1
                break  # only one mutation per read, so break after finding it

mutation_matrix.to_csv("MSQ0314 Plasmid/1_mutation_count_matrix.csv")
print(flag)