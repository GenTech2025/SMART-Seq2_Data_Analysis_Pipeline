# Identifying Gene Classes Exhibiting Stability in One Human Cell Type and Instability in Another to Unveil Biological Differences Between Cell Types
## Abstract
This project was a part of my master's degree and which revisits and extends upon a previous work conducted at the Tomlinson laboratory, where changes in the expression levels of mouse stably expressed genes (mSEGs) across different cell types in a mouse single-cell RNA sequencing dataset was used to uncover biological differences between two conditions. The first objective of this project was to re-evaluate and replicate the findings from the previous study, followed by an investigation into discrepancies in PCA results generated using Seurat’s RunPCA() function. Finally, we implemented a refined version of the previous algorithm on two human single-cell RNA sequencing datasets, processed in-house, to uncover the intricate transcriptional changes occurring between human embryonic stem cells and cells from early human development. The validity of the previous method was established through the successful replication of findings from the previous year and the PCA investigation revealed suceptibility of Seurat’s RunPCA function to implicit and explicit feature defination, compared to PCA functions from other packages. Finally, the implementation of the revised algorithm on human datasets revealed the intricate transcriptional changes that occur when embryonic stem cells transition from a pluripotent state to differentiated state. This analysis also highlighted how different genes exhibit stability at various stages of development, becoming unstable at other stages, thereby uncovering the complex biological differences between cell types.

<h2>Different Stages of the Project</h2>

<p align="center">
<img src="Final_Plots/dissertation-stage-overview.png" alt="Alt text" width="600">
</p>

<p>
  This project began as an exploratory study, with the primary goal of determining if an improved and more refined version of an algorithm previously applied to mouse datasets in a study at the <a href="https://regenerative-medicine.ed.ac.uk/research/simon-tomlinson">Tomlinson Laboratory</a> could be adapted for human single-cell RNA sequencing datasets. However, as the project progressed, new questions emerged, leading to its division into four major parts, each addressing a specific question and providing a foundation for the next steps.<br>
  <ol>
    <li>Part One: Re-evaluation and replication of findings from the previous study, along with investigation into discrepancies in PCA results generated using Seurat’s RunPCA() function.</li>
    <li>Part Two: Analysis of raw FASTQ files from a human single-cell RNA sequencing dataset generated using the SMART-Seq2 protocol.</li>
    <li>Part Three: Implementation of a refined version of the previous algorithm on two human single-cell RNA sequencing datasets. Analysis and interpretation of transcriptional changes between human embryonic stem cells and cells from early human development.</li>
  </ol>
</p>


<h2>Novel Algorithm To Find the Transcriptomics Differences between Two Cell Types/Groups/Stages</h2>
<p>
  The image below illustrates the different steps of the novel algorithm, which utilizes the concept of stably expressed genes as defined by <a href="https://doi.org/10.1093/gigascience/giz106">Yingxin et al 2019</a>. It combines gene set enrichment analysis to identify the intricate transcriptional changes that occur at the fundamental level. These changes shift the expression of the most stable genes between two conditions, providing insights into the key differences between the conditions.
</p>

<p align="center">
<img src="Final_Plots/fgsea-algorithm.png" alt="Alt text" width="600">
</p>

