# **Identification of Bottleneck Genes (Master Regulators) in Alzheimer's Disease Using Systems Biology Approaches**

## **Background**

Building on the work of [Dr. Andrea Califano’s lab (Columbia University)](http://califano.c2b2.columbia.edu/cancer-systems-biology), which focuses on the bottleneck hypothesis in cancer, this project seeks to apply similar concepts to Alzheimer’s disease. The bottleneck hypothesis indicates that although different genetic events may lead to a uniform disease phenotype, their effects converge on a small set of master regulator (MR) genes that control disease-driving processes within highly connected regulatory networks. These MRs, while not necessarily mutated, act as key bottlenecks, regulating the pathways and cellular programs that maintain disease states. Identifying these bottlenecks can provide mechanistic insights into disease biology and uncover potential therapeutic vulnerabilities.

<div style="text-align: center;">
  <img src="https://github.com/Gunja2502/Identifing-Bottleneck-genes-master-regulators-in-Alzheimer-s-Disease-RNA-seq-dataset/blob/e97a09a81ac85e2c792997ce0311951e31e46b73/images/Bottleneck%20genes%20(1).jpg" alt="Bottleneck genes (master-regulators)">
</div


Domain Selection
- Selected Domain: Regulatory Network Construction and Analysis
- Specific Application: Identification of Bottleneck genes (master regulators) associated with Alzheimer's Disease (AD)

---

## **Project Overview**

This project leverages systems biology tools to identify and rank key bottleneck genes, specifically transcription factors, acting as master regulators (MRs) in Alzheimer's Disease (AD). By applying **ARACNe (Algorithm for the Reconstruction of Accurate Cellular Networks)** and **VIPER (Virtual Inference of Protein Activity by Enriched Regulon analysis)**, we constructed gene regulatory networks and inferred protein activity levels from RNA-seq data.

---

## **Problem Statement**

Alzheimer's Disease is a progressive neurodegenerative disorder. Despite identifying numerous associated genes, understanding how these genes interact within regulatory networks remains a challenge. The goal of this project is to:

1. Construct and analyze gene regulatory networks to uncover regulatory bottlenecks.
2. Infer protein activity levels and identify potential therapeutic targets.
3. Identify transcription factors that serve as master regulators in AD.

---

## **Data Source**

- **Dataset:** [GSE153873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153873)  
  - **Title:** *An integrated multi-omics approach identifies epigenetic alterations associated with Alzheimer’s disease (RNA-seq with ERCC spike-in).*  
  - **Sample Groups:**  
    - **AD Group:** Patients with Alzheimer's Disease  
    - **Control Old Group:** Elderly individuals without AD  
    - **Control Young Group:** Younger individuals without AD  
  - **Gene Count:** 27,136 genes across 30 samples  

---

## **Tools and Technologies**

- **R Programming Language**
- **ARACNe-AP:** Network construction based on mutual information
- **VIPER:** Protein activity inference
- **Cytoscape:** Visualization and network analysis
- **Ensembl BioMart:** Gene length data extraction
- **TRRUST Database:** Transcription factor list

---

## **Project Workflow**
<img src="https://github.com/Gunja2502/Identifing-Bottleneck-genes-master-regulators-in-Alzheimer-s-Disease-RNA-seq-dataset/blob/main/images/Flowchart.jpg" alt="Flowchart" width="30%">

### 1. **Data Preprocessing**

- **Download Data:** Gene expression data in star counts format from GEO (GSE153873).
- **Log10 TPM Transformation:** Normalize the data by transforming it into TPM values.
- **Z-Score Standardization:** Standardize gene expression values for comparability.

### 2. **Gene Length Extraction**

- Extracted gene lengths from the **Ensembl BioMart** database to facilitate TPM calculation.

### 3. **Gene Regulatory Network Construction (ARACNe-AP)**

- **Steps:**
  1. Calculate mutual information (MI) between transcription factors and target genes.
  2. Apply Data Processing Inequality (DPI) to remove indirect interactions.
  3. Perform 100 bootstrap iterations for robustness --> Consolidate.
  4. Generate a final network file.

### 4. **Network Visualization (Cytoscape)**

- Visualize the ARACNe-constructed network.
- Perform **betweenness centrality** analysis to identify possible bottleneck genes.

### 5. **Protein Activity Inference (VIPER)**

- **Input:** Z-score transformed gene expression data and ARACNe network.
- **Output:** Protein activity scores for transcription factors.
- Visualize results as heatmaps for top regulators.

---

## **Key Files**

| **File**                       | **Description**                                               |
| ------------------------------- | ------------------------------------------------------------- |
| `GSE153873_summary_count.star.txt` | Raw gene expression data from GEO.                           |
| `mart_export.txt`              | Gene lengths extracted from Ensembl BioMart.                 |
| `log10_tpm_matrix.txt`         | Log10-transformed TPM gene expression data. Preprocessed gene expression data for ARACNe and VIPER.              |
| `network.txt`                  | ARACNe-AP output: inferred gene regulatory network.          |
| `Heatmap_Viper_Matrix2.pdf`    | Heatmap of inferred protein activity scores (VIPER output).  |
| `Heatmap_Gene_exp.pdf`         | Heatmap of gene expression data.                             |

---

## **How to Run the Project**

### 1. **Dependencies**

Install the required R libraries:

```r
install.packages(c("dplyr", "viper", "pheatmap", "org.Hs.eg.db", "ggplot2"))
```

## **2. Install Tools**

- **ARACNe-AP:** [ARACNe-AP GitHub Repository](https://github.com/califano-lab/ARACNe-AP)  
- **Cytoscape:** [Download Cytoscape](https://cytoscape.org/)

---

## **3. Run ARACNe-AP**

Execute the following command in PowerShell:

```bash
java -Xmx5G -jar aracne.jar -e "modified_file.txt" -o project --tfs "trrust_rawdata.human.txt" --pvalue 1E-8 --seed 1 --calculateThreshold
```

## **4. Visualize Network in Cytoscape**

- Import `network.txt` into Cytoscape.

---

## **5. Run VIPER Analysis**

- Run the VIPER R script to infer protein activity and generate heatmaps.

---

## **Results**

### **Visualised and analysed betweenness centrality of network generated by ARACNe-AP**
<div style="display: flex; justify-content: center; gap: 20px;">
  <img src="https://github.com/Gunja2502/Identifing-Bottleneck-genes-master-regulators-in-Alzheimer-s-Disease-RNA-seq-dataset/blob/main/images/Betweenness%20centrality%20.png" alt="Betweenness Centrality" width="50%">
  <img src="https://github.com/Gunja2502/Identifing-Bottleneck-genes-master-regulators-in-Alzheimer-s-Disease-RNA-seq-dataset/blob/main/images/Cytoscape%20network.png" alt="Cytoscape Network" width="41%">
</div>

Bottlenecks are proteins that have a high betweenness centrality, meaning they are network nodes that have many "shortest paths" going through them. Bottlenecks are key connector proteins that are more likely to be essential proteins. However, from the data/network we couldn't draw any conclusion as there were high number of genes in a single betweenness centrality value. 

### **Protein Activity Inference**

Regulators such as **CELSR3**, **CDAN1**, and **DMPK** show potential roles in AD.

### **Heatmaps**
<img src="https://github.com/Gunja2502/Identifing-Bottleneck-genes-master-regulators-in-Alzheimer-s-Disease-RNA-seq-dataset/blob/main/images/VIPER%20Heatmap.png" alt="Cytoscape Network" width="70%">

Visualizations of protein activity scores and gene expression patterns.

---

## **Future Scope**

### **Machine Learning Integration**

- Use Random Forest to identify and rank master regulators.

### **Experimental Validation**

- Validate candidate regulators in **in vitro** AD models.

### **Multi-Omics Integration**

- Incorporate proteomics and epigenomics data.

---

## **References**

1. **Califano Lab - Cancer Systems Biology:** [Link](http://califano.c2b2.columbia.edu/cancer-systems-biology)  
2. **ARACNe-AP GitHub:** [Link](https://github.com/califano-lab/ARACNe-AP)  
3. **VIPER Algorithm:** Alvarez et al., *Nature Genetics* (2016): [DOI](https://doi.org/10.1038/ng.3593)  
