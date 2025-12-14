# Plan to improve the project organization and presentation.

### Why make any changes to this project?
This project was a part of my dissertation for my MSc Degree in Data Science for Biology from The University of Edinburgh, when I first started to work on this project I had limited knowledge of good data science practices and was more focused on completing the project rather than focusing much on how the project was presented. Doing this was necessary at that time as time was limited and until the analysis was complete I was uncertain regarding the outcome of the project. However, not presenting it properly has prevented me from properly showing the effort and work put into this project and hence I have decided to reorganize and present the project in a more appropriate fashion.

### What changes to make?
- Re-organize the project directory structure, by focusing mainly on the implementation of scRNA seq pipeline and move the investigation regarding PCA discrepancy in the miscellaneous section.
- Rewrite the upstream stage of the pipeline (from raw FASTQ files to counts matrix) using Nextflow (hypothetical, dont have access to a HPC to verify its execution)
- Rewrite the downstream analysis script i.e. implementation of the novel algorithm in both R Markdown format and in a modularized script.
- Create comprehensive documentation and a updated readme with instructions on how to reproduce the results.
