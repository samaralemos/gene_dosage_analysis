# gene_dosage_analysis

This repository provides a Python-based pipeline for analyzing A-to-G gene dosage ratios in polyploid plant genomes. The workflow processes expression data from tetraploid and hexaploid species, calculates A/G ratios for diads and triads, performs statistical tests, and generates visualizations and summary reports.

The pipeline performs the following steps:

1. Calculates A/G gene expression ratios for diads (tetraploid) and triads (hexaploid) using normalized expression values (e.g., TPM).
2. Filters genes based on a minimum expression threshold of 0.5.
3. Compares A/G ratios between ploidy levels using the Kruskal-Wallis test and Dunn’s post hoc test.
4. Generates plots and statistical summaries.

Input Requirements

* `triads_file` (TSV): A file listing gene IDs for A and G subgenomes in tetraploids, and Am, At, and G subgenomes in hexaploids. The expected column order is:

  ```
  Am_gene_id    At_gene_id    G_gene_id
  ```

* `tetraploid_file` (TSV): Expression data file for the tetraploid genome, with gene IDs and corresponding expression values. 

* `hexaploid_file` (TSV): Expression data file for the hexaploid genome, including Am, At, and G gene IDs and their expression values.

How to Run:
To execute the analysis, edit the paths in the `main()` function near the bottom of the script:

```python
triads_file = "/path/to/triads_file.tsv"
tetraploid_file = "/path/to/tetraploid_expression.tsv"
hexaploid_file = "/path/to/hexaploid_expression.tsv"
output_folder = "gene_dosage_OUTPUTNAME"
main(triads_file, tetraploid_file, hexaploid_file, output_folder)
```

Estimated Runtime: For datasets containing around 10,000 triads, the complete analysis typically finishes in less than one minute on a standard desktop or laptop with 8 cores and 16 GB RAM. Larger datasets may require more time.

Output Files - The following files are saved to the specified output folder:

| File Name          | Description                                                |
| ------------------ | ---------------------------------------------------------- |
| diads\_ratios.tsv  | A/G ratios computed for tetraploid diads                   |
| triads\_ratios.tsv | A/G ratios computed for hexaploid triads                   |
| intersect.tsv      | Merged data for matching diad and triad entries            |
| statistics.txt     | Results of Kruskal-Wallis and Dunn's post hoc tests        |
| summary.txt        | Mean, median, and standard deviation of A/G ratios         |
| plots.pdf          | Boxplots of A/G ratios and normalized triad-to-diad ratios |

