from collections import Counter

reference_seq = "MPHSSLHPSIPCPRGHGAQKAALVLLSACLVTLWGLGEPPEHTLRYLVLHLASLQLGLLLNGVCSLAEELRHIHSRYRGSYWRTVRACLGCPLRRGALLLLSIYFYYSLPNAVGPPFTWMLALLGLSQALNILLGLKGLAPAEISAVCEKGNFNV"#seg1
#reference_seq = "AHGLAWSYYIGYLRLILPELQARIRTYNQHYNNLLRGAVSQRLYILLPLDCGVPDNLSMADPNIRFLDKLPQQTGDRAGIKDRVYSNSIYELLENGQRAGTCVLEYATPLQTLFA" #seg2
#reference_seq = "MSQYSQAGFSREDRLEQAKLFCRTLEDILADAPESQNNCRLIAYQEPADDSSFSLSQEVLRHLRQEEKEEVTVGSLKTSAVPSTSTMSQEPELLISGMEKPLPLRTDFS" #seg3

input_file = "MSQ0314 Plasmid/plas1_translated.txt"
output_file = "MSQ0314 Plasmid/1_one_mut.txt"

mutation_counts = Counter()

with open(input_file, "r") as f, open(output_file, "w") as fout:
    for line_num, line in enumerate(f, 1):
        read_seq = line.strip()
        if len(read_seq) != len(reference_seq):
            # Skip sequences with incorrect length
            print(f"Warning: line {line_num} length {len(read_seq)} != reference length {len(reference_seq)}")
            continue

        # Count mutations (differences)
        mutations = sum(1 for a, b in zip(read_seq, reference_seq) if a != b)
        mutation_counts[mutations] += 1

        if mutations == 1:
            fout.write(read_seq + "\n")

with open("MSQ0314 Plasmid/1_counts.txt", "w") as f:
    for mut_count in sorted(mutation_counts.keys()):
        f.write(f"{mutation_counts[mut_count]}\n")

print("Mutations\tRead_count")
for mut_count in sorted(mutation_counts.keys()):
    print(f"{mut_count}\t{mutation_counts[mut_count]}")

ss = 0
for i in range(len(mutation_counts)):
    ss += mutation_counts[i]
print(ss)