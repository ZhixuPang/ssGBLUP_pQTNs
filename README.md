# ssGBLUP_pQTNs

Sure! Here's an English version of the README file for your GitHub project:

---

# Genomic Prediction Analysis Example

This repository provides an example of genomic prediction analysis using genotype and phenotype data. The analysis includes calculations of genetic relationship matrices (A, G) and performs a mixed model analysis to predict breeding values.

The code demonstrates how to read PLINK format data, calculate the genomic relationship matrices, perform genome-wide association studies (GWAS), select pQTNs (pseudo Quantitative Trait Nucleotides), and finally use a mixed model to predict breeding values.

## Installation

To run this code, you need to install the following Julia packages:

```julia
using Pkg
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Random")
Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
Pkg.add("Plots")
Pkg.add("StatsBase")
```

## Data Requirements

This analysis requires several data files in specific formats:

1. **Genotype Data**: PLINK format genotype data, which includes `.fam`, `.bim` and `.fam` files.
2. **Phenotype Data**: A text file (`phenotypes.txt`) containing individual IDs and phenotype values.
3. **Pedigree Data**: A pedigree file (`ped.txt`) containing family relationships.

Example file formats:

**PLINK Files (`qtl-mas-2010.ped` and `qtl-mas-2010.map`)**:
- The `.ped` file contains the genotype data for individuals in the following format:
    ```
    FID IID FatherID MotherID Sex Phenotype
    1   1   0       0       1   2
    1   2   0       0       2   3
    ```

**Phenotype Data (`phenotypes.txt`)**:
- The phenotype file contains individual IDs and their corresponding trait values:
    ```
    ID  Q   B
    1   1.2 2.3
    2   3.4 5.6
    ```

**Pedigree Data (`ped.txt`)**:
- The pedigree file contains family relationships for each individual:
    ```
    FID IID FatherID MotherID Sex
    1   1   0       0       1
    1   2   0       0       2
    ```

## Workflow

### Step 1: Read Genotype Data

The genotype data is read from PLINK files using the `readPlink` function. The `geno` object contains the genotype matrix along with individual and marker information.

```julia
geno = readPlink("qtl-mas-2010")
```

### Step 2: Random Sampling of Individuals

2000 individuals are randomly selected from the genotype data for further analysis:

```julia
Random.seed!(1234)
gid = sort(sample(geno.fam.Individual_ID, 2000, replace=false))
```

### Step 3: Calculate Allele Frequencies and Genotype Data

The allele frequencies are calculated, and individuals with low allele frequencies are excluded:

```julia
p = mean(geno.geno, dims=1)' ./ 2
p = p[:, 1]
geno = Geno(geno.geno[gid, findall(p .> 0.05)], copy(geno.fam[gid, :]), copy(geno.map[findall(p .> 0.05), :]))
```

### Step 4: Calculate Pedigree and Genetic Relationship Matrices

The A matrix (pedigree relationship matrix) and the G matrix (genomic relationship matrix) are computed:

```julia
ped = getPedigree("ped.txt"; separator = " ")
A = calA(ped)
G = calG(geno)
```

### Step 5: Perform GWAS and Select pQTNs

A genome-wide association study (GWAS) is performed, and pQTNs are selected for inclusion in the prediction model:

```julia
meta = gwas(pheno.Q[gwasid], 
            Geno(geno.geno[1:length(gwasid), :], copy(geno.fam[1:length(gwasid), :]), copy(geno.map[:, :])), 
            hcat(ones(length(gwasid), 1), pca10), 
            Symmetric(G[1:length(gwasid), 1:length(gwasid)]), var(pheno.Q) * 0.5, 1.)
```

### Step 6: Mixed Model Prediction

Using the selected pQTNs, a mixed model is solved to predict breeding values for the individuals:

```julia
b, u = solveMME(X[pheno.seq_IDs, :], Z[pheno.seq_IDs, :], Hinv, 1., pheno.Q)
yhat = X * b + u
```

### Step 7: Evaluation

The predicted breeding values (`yhat`) are compared to the true breeding values (`tbv`) by calculating the correlation:

```julia
cor(yhat[2327:end], tbv.Q_Tbv[2327:end])
```

## Results

The code predicts breeding values based on genomic and pedigree information. The final output includes predicted values and the correlation with true breeding values.

## Notes

- Ensure all input files are in the correct format.
- Modify parameters (e.g., weights, regularization) to optimize prediction accuracy.
- GWAS (Genome-Wide Association Studies) can be performed using various software tools like GEMMA and GCTA.
