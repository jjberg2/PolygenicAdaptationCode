CreateTraitFile <- function ( gwas.data.filename , gen.data ) {

	recover()	
	# trait.effects.data is a data table that was already read in by R in the EstimateTraitEffects function call. In that function, the name of a file is given, and it is read in as trait.effects.data. In the CreatTraitFiles call, it is given the same name, and passed in as a data table.
	# gen.data is passed in as a string, the name of a file, from the EstimateTraitEffect call. The file it specifies should have columns of SNP, CLST, FRQ, A1, A2
	# strand.data is a string specifying the name of a file containing columns of SNP, and STRAND (rsID, and "+" or "-")
	# output.name specifies the name of the file to be output containing the extracted SNP frequencies
	# temp <- unlist ( strsplit ( gen.data , split = "/" ) )
	# gen.data.name <- temp [ length ( temp ) ]
	gwas.data <- read.table ( gwas.data.filename , header = T , stringsAsFactors = F )
	temp <- unlist ( strsplit ( gen.data , split = "/" ) )
	temp <- temp [ length ( temp ) ]
	
	write ( paste ( gwas.data$SNP , sep = "\n" ) , file = paste ( gwas.data.filename , ".search" , sep = "" ) , ncolumns = 1 )
	my.cmd <- paste ( "Scripts/sampleSNPs.pl ", gwas.data.filename , ".search " , gen.data , " > ", gwas.data.filename , "." , temp , ".freqs.temp" , sep = "")
	system ( my.cmd )

	# The gen.data file I'm currently using has a number of entries with asterisks after the base ID, which just indicate that the derived allele was not known when Joe Pickrel did the imputation.
	system ( paste ( "sed 's/*//'" , " " , gwas.data.filename , "." , temp , ".freqs.temp" , " > " , gwas.data.filename , "." , temp , ".freqs" , sep = "" ) )
	system ( paste ( "rm -f " , gwas.data.filename , "." , temp , ".freqs.temp" , sep = "" ) )	
	
	
}