# MIDESP

Mutual information based detection of epistatic SNP pairs

MIDESP calculates the mutual information between SNP pairs and phenotypes and reports the top interactions.

Genotype and phenotype data needs to be in the `tped` and `tfam` format used by [PLINK](https://www.cog-genomics.org/plink/1.9/formats).

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
## Examples
The datasets we analyzed in our publication can be found under [Tuberculosis](https://doi.org/10.5061/dryad.519bm) and [Eggweight](https://figshare.com/articles/dataset/Genome-wide_Association_Analysis_of_Age-Dependent_Egg_Weights_in_Chickens/5844420) (see also their corresponding publications [https://doi.org/10.1038/hdy.2013.137](https://doi.org/10.1038/hdy.2013.137) and [https://doi.org/10.3389/fgene.2018.00128](https://doi.org/10.3389/fgene.2018.00128)).

##### Data preparation 

After converting them to a PLINK readable format we first filtered them and performed SNP pruning with PLINK using the following commands.

**Eggweight**

This dataset contains multiple phenotypes. We have only analyzed the phenotype **EW36**.

```
plink --allow-no-sex --geno 0.03 --hwe 0.000001 --maf 0.01 --mind 0.05 --prune --recode transpose --chr-set 35 --out Eggweight_Filtered --tfile Eggweight
plink --allow-no-sex --indep-pairwise 10000 5 0.99 --make-founders --chr-set 35 --out Eggweight_Filtered_PruneInfo --tfile Eggweight_Filtered
plink --allow-no-sex --extract Eggweight_Filtered_PruneInfo.prune.in --make-founders --recode transpose --chr-set 35 --out Eggweight_Filtered_Pruned --tfile Eggweight_Filtered
```

**Tuberculosis**

This dataset contains multiple SNPs with the same ID which we remove in the first step (IDs are given in the file Tuberculosis.dupIDs in TuberculoseData.tar.gz) as well as several SNPs with chr 0, which we need to add again after SNP pruning.

```
plink --allow-no-sex --geno 0.03 --hwe 0.000001 --maf 0.01 --mind 0.05 --exclude Tuberculosis.dupIDs --prune --recode transpose --cow --out Tuberculosis_Filtered --tfile Tuberculosis
plink --allow-no-sex --chr 0 --recode --cow --out Tuberculosis_Filtered_Chr0 --tfile Tuberculosis_Filtered
plink --allow-no-sex --indep-pairwise 10000 5 0.99 --cow --out Tuberculosis_Filtered_PruneInfo --tfile Tuberculosis_Filtered
cut -f2 Tuberculosis_Filtered_Chr0.map >> Tuberculosis_Filtered_PruneInfo.prune.in
plink --allow-no-sex --extract Tuberculosis_Filtered_PruneInfo.prune.in --recode transpose --cow --out Tuberculosis_Filtered_Pruned --tfile Tuberculosis_Filtered

```

Finally, for storing them on Github we converted them to the binary format used by PLINK.

```
plink --allow-no-sex --chr-set 35 --out Eggweight_Filtered_Pruned --make-bed --tfile Eggweight_Filtered_Pruned
plink --allow-no-sex --cow --out Tuberculosis_Filtered_Pruned --make-bed --tfile Tuberculosis_Filtered_Pruned
```

These final datasets can be found in EggweightData.tar.gz and TuberculoseData.tar.gz. 

They first need to be converted to the `tped` and `tfam` format before MIDESP can analyze them.

```
plink --allow-no-sex --chr-set 35 --out Eggweight_Filtered_Pruned --recode transpose --bfile Eggweight_Filtered_Pruned
plink --allow-no-sex --cow --out Tuberculosis_Filtered_Pruned --recode transpose --bfile Tuberculosis_Filtered_Pruned
```

##### Analysis with MIDESP

In our analysis we then used the following commands:

```
java -jar MIDESP.jar -threads 70 -out Eggweight_Filtered_Pruned.epi -keep 0.25 -fdr 0.005 -apc 5000 -cont -k 30 Eggweight_Filtered_Pruned.tped Eggweight_Filtered_Pruned.tfam
java -jar MIDESP.jar -threads 70 -out Tuberculosis_Filtered_Pruned.epi -keep 0.1 -fdr 0.005 -apc 5000 Tuberculosis_Filtered_Pruned.tped Tuberculosis_Filtered_Pruned.tfam 
```

## License

This project is licensed under the **GPL-3.0 License** - see [LICENSE](LICENSE) for more information.