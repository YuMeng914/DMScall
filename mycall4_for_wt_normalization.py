import csv
import numpy as np
import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap

input_file = "MSQ0333/1l_mutation_count_matrix.csv"
ctrl = []
with open(input_file, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader:
        values = row[1:21]
        int_values = [int(float(val)) if val.strip() else 0 for val in values]
        ctrl.extend(int_values)

input_file = "MSQ0333/1h_mutation_count_matrix.csv"
high = []
with open(input_file, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader:
        values = row[1:21]
        int_values = [int(float(val)) if val.strip() else 0 for val in values]
        high.extend(int_values)

ratio = []
sc = sum(ctrl)
sh = sum(high)

with open("MSQ0333/1l_counts.txt", "r") as f:
    low_readcounts = [float(line.strip()) for line in f]
with open("MSQ0333/1h_counts.txt", "r") as f:
    high_readcounts = [float(line.strip()) for line in f]
norm_wt = math.log2((high_readcounts[0]/sum(high_readcounts)) / (low_readcounts[0]/sum(low_readcounts)))
print(norm_wt)

threshold = 20
for i in range(len(ctrl)):
    if ctrl[i] > threshold and high[i] > threshold:
        ratio.append(math.log2((high[i]/sh) / (ctrl[i]/sc)) - norm_wt)
    else:
        #ratio.append(0)
        ratio.append(np.nan)

wt_sequence = "MPHSSLHPSIPCPRGHGAQKAALVLLSACLVTLWGLGEPPEHTLRYLVLHLASLQLGLLLNGVCSLAEELRHIHSRYRGSYWRTVRACLGCPLRRGALLLLSIYFYYSLPNAVGPPFTWMLALLGLSQALNILLGLKGLAPAEISAVCEKGNFNV"#seg1
#wt_sequence = "AHGLAWSYYIGYLRLILPELQARIRTYNQHYNNLLRGAVSQRLYILLPLDCGVPDNLSMADPNIRFLDKLPQQTGDRAGIKDRVYSNSIYELLENGQRAGTCVLEYATPLQTLFA" #seg2
#wt_sequence = "MSQYSQAGFSREDRLEQAKLFCRTLEDILADAPESQNNCRLIAYQEPADDSSFSLSQEVLRHLRQEEKEEVTVGSLKTSAVPSTSTMSQEPELLISGMEKPLPLRTDFS" #seg3
aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

ratio = np.array(ratio).reshape((len(wt_sequence),20))

dff = pd.DataFrame(ratio)

wt_mask = np.zeros_like(ratio, dtype=bool)
for i, wt_aa in enumerate(wt_sequence):
    if wt_aa in aa_list:
        j = aa_list.index(wt_aa)
        wt_mask[i, j] = True

nan_mask = np.isnan(ratio)

ratio = ratio.T
nan_mask = nan_mask.T
wt_mask = wt_mask.T

def truncate_colormap(cmap, minval=0.1, maxval=0.9, n=256):
    new_cmap = LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

original_cmap = plt.get_cmap("RdBu_r")
lighter_cmap = truncate_colormap(original_cmap, 0.1, 0.9)
lighter_cmap.set_bad(color='#D3D3D3')

plt.figure(figsize=(30,4))
#cmap = sns.color_palette("RdBu_r", as_cmap=True)
#cmap.set_bad(color='#D3D3D3')
hmp = sns.heatmap(ratio, square=True, mask=nan_mask, annot=False, cmap=lighter_cmap, linewidths=0.25, linecolor='gray', vmin=-1, vmax=1, cbar_kws={"shrink": 1}, yticklabels=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

for i in range(20):
    for j in range(len(ratio[0])):
        if wt_mask[i, j]:
            plt.plot(j + 0.5, i + 0.5, 'o', color='#414A4C', markersize=4, zorder=2)

tick_pos = [0,9,19,29,39,49,59,69,79,89,99,109,119,129,139,149]
tick_lab = [1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
#Section1

#tick_pos = [4,14,24,34,44,54,64,74,84,94,104,114]
#tick_lab = [160,170,180,190,200,210,220,230,240,250,260,270]
#Section2

#tick_pos = [9,19,29,39,49,59,69,79,89,99]
#tick_lab = [280,290,300,310,320,330,340,350,360,370]
#Section3

#tick_pos = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375]
#tick_lab = [1, 26, 51, 76, 101, 126, 151, 176, 201, 226, 251, 276, 301, 326, 351, 376]
#FL

tick_pos = [i + 0.5 for i in tick_pos]
hmp.set_xticks(tick_pos)
hmp.set_xticklabels(tick_lab)
savepath = 'MSQ0333/seg1_normWT_thre20_final.pdf'
plt.tick_params(axis='y', length=0)
fon = fm.FontProperties(family="Arial")
plt.yticks(rotation=0, va='center', ha='center', fontproperties=fon, fontweight='bold')
plt.xticks(rotation=0, ha='center', fontproperties=fon, fontweight='bold')
plt.tight_layout()
plt.savefig(savepath, format='pdf', dpi=300)
plt.show()
