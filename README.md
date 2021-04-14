# AMscripts_mouseEmbryoFISH

GitHub folder containing scripts to 
1) Generate data for imputation for SpatialEmbryos and analyze the performance.
2) Part of MHB analysis (virtual dissection of MHB and differential expression between Midbrain and Hindbrain regions).
3) Get intensity thresholds for smFISH genes and get counts matrices for smFISh genes.

`Scripts:`

`.\core_functions.R\` contains functions that are frequently used throughout other pipleines.

`.\visualization_functions.R\` contains color codes for celltypes.

---

**1. Imputation for Spatial embryos:**

a. `.\generateData\imputation\mapping\` contains scripts to get mapping between seqFISH and scRNA-seq; scRNA-seq onto itself.

b. `.\generateData\imputation\performance\` contains scripts to get intermediate imputations (for each seqFISH gene) and accordingly intermediate prediction scores.

c. `.\generateData\imputation\imputation\` contains scripts to get final imputations and prediction scores for each gene / embryo / z-slice.

d. `.\analysis\perfromance_intermediate_imputation\` contains script to analyse the intermediate imputations and comapre them against experimentally measured seqFISH data.

e. `.\analysis\imputation_gene_prediction\` contains script to combine prediction scores for each gene / embryo/ z-slice for the Supplementary table.

f.  `.\analysis\comparison_w_smFISH\` contains script to compare final imputations across independent experimental validation (smFISH).

---

**2. Midbrain-Hindbrain border analysis:**

a. `.\analysis\MHB\` contains script to virtually dissect MHB region and perform differential expression analysis.

---

**3. smFISH data:**

a. `.\analysis\smFISH_channel_effect\` contains script to predict intensity thresholds for each gene (that was probed for smFISH)/ field of view.

b. Based on predicted intensity thresholds, the counts matrix was retrieved using custom MATLAB script (in the same manner that for genes probed for the main experiment, availbale upon request).

c. `.\generateData\sce_smFISH_genes\` contains script to generate SingleCellExperiment object from counts matrix.

---
