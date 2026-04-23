# OSCA Basics - Chapter 2: Normalization (Comprehensive Study Notes)

## 1. Mathematical Notation & Glossary
Before diving into formulas, we define the variables used in this chapter:

* **$i$**: The index for a **Cell** (the column in your data matrix).
* **$g$**: The index for a **Gene** (the row in your data matrix).
* **$C_{ig}$**: The **Raw Count** (number of reads or UMIs) for gene $g$ in cell $i$. This is the value found at the intersection of a specific row and column.
* **$L_i$**: The **Library Size** of cell $i$. It is the sum of all raw counts for that cell ($\sum_{g} C_{ig}$).
* **$N$**: The **Total Number of Cells** in the entire dataset.
* **$s_i$**: The **Size Factor** for cell $i$. This is the scaling constant we calculate to remove bias.
* **$X_{ig}$**: The **Normalized Expression Value** for gene $g$ in cell $i$ after scaling.
* **$logcounts_{ig}$**: The final **Log-transformed** normalized value.
* **Pseudocount**: A small constant (usually **+1**) added to counts before log-transformation to avoid calculating the log of zero.

---

## 2. The Core Motivation: Why Normalize?
Systematic technical differences (cDNA capture, PCR efficiency, sequencing depth) create "noise" that varies cell-to-cell.
* **Technical Bias:** A cell might have higher counts simply because it was sequenced more deeply.
* **Goal:** Ensure that differences in expression profiles are driven by **Biology**, not technical biases.

---

## 3. The Size Factor ($s_i$)
The Size Factor is the "scaling ruler" for each cell. It represents the relative technical bias of that cell compared to the rest of the population.

### How to Calculate a Library Size Factor ($s_i$):
For a specific cell $i$, the size factor is its library size divided by the average library size of all cells:

$$s_i = \frac{L_i}{\frac{1}{N} \sum_{j=1}^{N} L_j}$$

* **Mean Centering:** This ensures the average of all $s_i$ in the dataset is exactly **1**.

### Applying the Scaling Formula:
To get the normalized value ($X_{ig}$), we divide the raw count by the size factor:
$$X_{ig} = \frac{C_{ig}}{s_i}$$

---

## 4. Worked Example: Calculating $s_i$
Dataset: **3 cells**
* **Cell 1 ($L_1$):** 2,000 counts
* **Cell 2 ($L_2$):** 5,000 counts
* **Cell 3 ($L_3$):** 8,000 counts

1.  **Calculate Mean Library Size:** $(2000 + 5000 + 8000) / 3 = \mathbf{5,000}$
2.  **Calculate $s_1$:** $2,000 / 5,000 = \mathbf{0.4}$
3.  **Calculate $s_2$:** $5,000 / 5,000 = \mathbf{1.0}$
4.  **Calculate $s_3$:** $8,000 / 5,000 = \mathbf{1.6}$

**Interpretation:** If Gene A has **16 counts** in Cell 3 ($C_{3,A} = 16$), its normalized value is $16 / 1.6 = \mathbf{10}$.

---

## 5. Scaling Normalization Methods

### A. Library Size Normalization
* Uses the total count ($L_i$) to find $s_i$.
* **Pitfall: Composition Bias.** If a few genes are extremely active, they "hog" the sequencing space, making other genes look downregulated even if they aren't.

### B. Normalization by Deconvolution (`scran`)
* **Method:** Groups cells into clusters and **pools** their counts to increase signal and reduce zeros.
* **Benefit:** Much more robust against composition bias; doesn't assume every cell has the same total mRNA content.

### C. Spike-in Normalization (ERCC)
* Uses synthetic "alien" RNA added at a fixed concentration.
* **Purpose:** To preserve biological changes in **Total RNA Content** (e.g., if a cell physically grows and triples its mRNA).

---

## 6. Log-Transformation: $log_2(x + 1)$
Final values are stored in the `logcounts` assay. 

### The Combined Formula:
$$\text{logcounts}_{ig} = \log_2\left(\frac{C_{ig}}{s_i} + 1\right)$$

### Why Log-Transform?
1.  **Focus on Fold-Changes:** $log_2$ makes a "doubling" equal to a distance of **1**.
2.  **Variance Stabilization:** It prevents high-abundance genes from dominating the analysis purely because of their large raw numbers.
3.  **Distance-based analysis:** Necessary for PCA, Clustering, and UMAP.

---

## 7. Bioconductor Implementation
* `librarySizeFactors(sce)` -> Simple scaling.
* `quickCluster(sce)` + `calculateSumFactors(sce)` -> Deconvolution.
* `logNormCounts(sce)` -> Applies scaling + log-transformation; creates `logcounts`.