# MIDESP

Mutual information based detection of epistatic SNP pairs

MIDESP calculates the mutual information between SNP pairs and phenotypes and reports the top interactions.

Genotype and phenotype data needs to be in the `tped` and `tfam` format used by [PLINK](https://www.cog-genomics.org/plink2/formats).

The program will create two files:

- `outputfile.sigSNPs` which contains a list of the SNPs that were found as being strongly associated to the phenotype along with their entropy and the mutual information between the SNP and the phenotype
- `outputfile.epi` which contains a list of the top SNP pairs that show the strongest association to the phenotype along with the mutual information between the SNP pair and the phenotype and the mutual information corrected through application of the APC theorem

## Requirements

Java 8 or later

[Apache Commons Math](https://github.com/apache/commons-math) 3.6.1 or later (for compiling the program yourself)

## Getting Started

1. Download MIDESP.jar from [Releases](https://github.com/FelixHeinrich/MIDESP/releases/tag/1.0) or compile it yourself from the source code
2. Open a terminal and go to the directory where MIDESP.jar is located
3. Run the following command with suitable parameters


```
java -jar MIDESP.jar {Options} tpedFile_Path tfamFile_Path
```

Optional parameters:

```
-out		file	name of outputfile (default tpedFile.epi)
-threads	number	number of threads to use (default = Number_of_Cores / 2)
-keep		number	keep only the top X percentage pairs with highest MI (default = 1)
-cont			indicate that the phenotype is continuous
-k		number	set the value of k for MI estimation for continuous phenotypes (default = 30)
-fdr		number	set the value of the false discovery rate for finding significantly associated SNPs (default = 0.005)
-apc		number	set the number of samples that should be used to estimate the average effects of the SNPs (default = 5000)
```

## License

This project is licensed under the **GPL-3.0 License** - see [LICENSE](LICENSE) for more information.