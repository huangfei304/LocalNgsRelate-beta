# LocalNgsRelate-beta

## Brief description

​        LocalNgsRelate-beta is based on [LocalNgsRelate]( https://github.com/idamoltke/LocalNgsRelate) and make some change for reducing memory, time and output format for following machine-learning (ML) analysis. Using the LocalNgsRelate-beta result and machine-learning model, to determine the kinship. 

​       LocalNgsRelate can be used to infer IBD sharing along the genomes of two individuals from low-depth Next Generation Sequencing (NGS) data by using genotype likelihoods (instead of called genotypes).  To be able to infer the IBD sharing you will need to know the population allele frequencies and have genotype likelihoods. This can be obtained by some software ( e.g. the program [ANGSD](http://popgen.dk/angsd/index.php/Quick_Start)).

## How to download and install
On a linux or mac system with git and g++ installed LocalNgsRelate can be downloaded and installed as follows:

```
1. LocalNgsRelate
   git clone https://github.com/idamoltke/LocalNgsRelate
   cd LocalNgsRelate/src/cpp/
   make
2. LocalNgsRelate-beta
   https://github.com/huangfei304/LocalNgsRelate-beta.git
   cd LocalNgsRelate-beta/src
   make
```

## Input file format
### 1. Formal description
LocalNgsRelate-beta takes three files as input: two files with genotype likelihoods (-a and -b) and a file with population allele frequencies (-f) for the sites there are genotype likelihoods for with the same index ( or order).

The genotype likelihood file needs to contain a line for each site with 3 values and it needs to be in beagle format and gz compressed (see e.g. http://www.popgen.dk/angsd/index.php/Beagle_input). Note that the marker name needs to be of the allele frequencies. The reason for this is that the programs needs to know which markes are on the same chromosome and the position of each marker.

The frequency file needs to contain a line per site with the allele frequency of the site in it. 
For examples of the two types of input files see the files in the folder exampledata which are described below.


### Example input files
The example files included here (in the folder test). There are two files:

1. test.beagle  a file with genotype likelihoods with 21 sites

   | marker | allele1 | allele2 | Ind0   | Ind0   | Ind0   |
   | ------ | ------- | ------- | ------ | :----- | ------ |
   | 100001 | 2       | 0       | 0.9846 | 0.0154 | 0.0000 |
   | 100002 | 3       | 1       | 0.0000 | 1.0000 | 0.0000 |
   | 100005 | 0       | 2       | 0.0000 | 0.0019 | 0.9981 |
   | 100006 | 0       | 2       | 0.9697 | 0.0303 | 0.0000 |

2. test.freq.gz       a file with allele frequency estimates for 1000 sites in han population (no header).

   | 0    | chr1_13273 | G    | C    | 0.111726 |
   | ---- | ---------- | ---- | ---- | :------- |
   | 1    | chr1_13649 | G    | C    | 0.113274 |
   | 2    | chr1_16288 | C    | G    | 0.285509 |
   | 3    | chr1_17020 | G    | A    | 0.14469  |
   | 4    | chr1_20144 | G    | A    | 0.132412 |



As an example of a pair of unrelated individuals you can e.g. use NA19027 and NA19313, so the samples with index 0 and 2. For a description of how this dataset was made see "Making input data" below.

## Output format
Successfully running the program should lead to 1 output file and if run the program with "–o test" these will be called

1) test.IBDmerge.gz

```
total	ibd1	ibd1R	ibd2	ibd2R	kinship	ibd1L	ibd2L
2808418773	759403950	0.06760066886	429275341	0.076426518923	0.1440271877	9752625	1934761

```


## Run examples

### Making input data
For getting the input data, just see the [LocalNgsRelate](https://github.com/idamoltke/LocalNgsRelate), which giving the detailed information for how the the beagle format individual result.

###	Analysing pair samples
To determine the IBD and related result. Creating a folder firstly and moving into this folder. Then run the following commands to analyse the pair of samples.

```
mkdir test
cd test
../localngsrelate-beta -a A.beagle.gz -b B.beagle.gz -f allchrs.freq.gz -o A_B
```

###	Compare LocalNgsRelate and LocalNgsRelate-beta

​      using all the common variant (maf >5%) in 1000 Genomes Project (1kGP), the memory and time used is as follows:  

| list                | memory | run time |
| ------------------- | ------ | -------- |
| LocalNgsRelate      | ~30G   | ~30min   |
| LocalNgsRelate-beta | ~4G    | ~2min    |

### Infer the degree of relatedness with ML

1) using the kinship coefficient to infer the degrees of relatedness for 1000 genomes project samples.	

   <img src="/ML/1kg_kinship.png" alt="1kg_kinship" style="zoom:50%;" />

2) using the kinship coefficient to infer the degrees of relatedness  for test samples.	

   <img src="/ML/lowpass_kinship.png" alt="lowpass_kinship" style="zoom:40%;" />

   Note: **PC**(parent-child) , **MZ**(monozygotic  twins),  **FS**(full sibling), **GG**(grandparent-grandchild), **AUNN**(Aunt-Uncle-Niece-Nephew), **HS**(half sibling), **1stC**(full first  cousin),  **UN** (>=2ndC). 




## Citing and references
1. Severson, A.L., Korneliussen, T.S. and Moltke, I., 2022. LocalNgsRelate: a software tool for inferring IBD sharing along the genome between pairs of individuals from low-depth NGS data. *Bioinformatics*, *38*(4), pp.1159-1161.

