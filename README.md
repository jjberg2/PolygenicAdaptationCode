PolygenicAdaptationCode
=======================
This repository contains R code implementing the methods described in my paper with Graham Coop, titled "A Population Genetic Signal of Polygenic Adaptation": 10.1371/journal.pgen.1004412


Here is some minimal documentation describing each of the files required to run the analysis. See the example run file at path/Example/Run_Files/exampleHeightRunFile.R for an example of a function call. All file names should be specified either with a full path or a path relative to the directory from which the function is called.

See also the headers of each of the corresponding example files found in the Example/ directory for a few more lines of documentation for each file type

gwas.data.file contains the effect size information for the SNPs identified in the GWAS

freqs.file contains the frequency information from the test panel for the SNPs ID'd in the GWAS

env.var.data.files contain the information on the environmental variables as well as groupings for the regional Z scores. Regional groupings should be given as integers. At least one of these files must be included, even if you are only interested in the general overdispersion test. In this case you can just fill in some arbitrary integers for the ENV and REG columns and then ignore those results

match.pop.file contains the appropriate match variables (frequency/imputation status/B Value/whatever other variable you may find useful) for all SNPs which may be included in the covariance matrix or for simulating null genetic values.

full.data.set.file file containing all of the frequency information for all populations

path is the path to the folder you want to write output to, relative to the directory from which you are running the analysis.

match.categories gives the column names for the variables you want to match when sampling null SNPs. MAF (minor allele frequency) or FRQ (if alleles are polarized by ancestry and you want to use this information in matching) is required, others are optional. 

match.bins gives instructions for how to construct the matching bins. For example, I split minor allele frequency into bins of 0.02 by giving the argument "seq ( 0 , 0.5 , 0.02 )", and I split the categorical variable for imputation status (either 0 or 1) by giving the argument c ( 2 ) to indicate that there were two bins.

cov.SNPs.per.cycle: R is memory limited, so when the construction of the covariance matrix can be split into multiple batches. In this case, I read 5000 SNPs per batch. The number of SNPs which can successfully be read is a function of the size of the dataset (# of populations), so you can just trial and error it, but my MacBook Air 11" with 4 GB RAM easily handled 5000; one could also calculate a covariance matrix with some other piece of software that's less memory limited than my script, and then just save it in the output folder under the name "uncenteredcovmat.Rdata" and set "load.cov.matrix" (see below) to TRUE so that it will be read in rather than calculated from my script. One could also calculate the covariance matrix using some other piece of software that wasn't as memory limited as my scripts, then save it as "uncenteredcovmat.Rdata" in the appropriate directory and set the "load.cov.matrix" flag to TRUE (see below) to just have my program load it in.

cov.cycles: Gives the number of batches to average the covariance matrix over

null.phenos.per.cycle: essentially the same as cov.SNPs.per.cycle, but specifies the number of null genetic values to be calculated rather than the number of SNPs to be sampled. For the height data, which included 160 SNPs for 52 populations for each set of null genetic values, my computer freaked out at anything other 1000 per batch, but for skin pigmentation, which had only 4 SNPs per genetic value, 10000 was a piece of cake.

null.cycles: number of batches of null genetic values to calculated to get the null distributions.

load.cov.matrix: if set to TRUE, the program will check in the folder path/Output for a file called "uncenteredcovmat.Rdata", which may have been generated in a previous run of the program ('path' here represents the path given as a part of the function call). If the file is present, it is read in, and calculation of the covariance matrix is skipped. If absent, or if the value of FALSE, the covariance matrix is calculated as normal. 

sim.null: if FALSE, the an analytical p-value is obtained for the Qx statistic only by comparing to the relevant Chi^2 distribution, and p values for the Z test are taken from the standard normal distribution, rather than from empirical resampling.

check.allele.orientation: compare the gwas.data.file and the freq.file to make sure that alleles are oriented the same way around by comparing columns A1 and A2 from the two files. Alleles which are polarized in opposite directions will be flipped (the frequency listed in the files is assumed in all cases to be that of the A1 allele). The program throws an error if a suspected strand issue is encountered. It does not check for strand errors if this flag is set to FALSE.



Copyright (C) 2014 by Jeremy J Berg <jeremy.jackson.berg@gmail.com>

Permission to use, copy, modify, and/or distribute this software for any purpose
with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
