CreateTraitFile <- function ( gwas.data.filename , gen.data ) {

	#recover()	
	# trait.effects.data is a data table that was already read in by R in the EstimateTraitEffects function call. In that function, the name of a file is given, and it is read in as trait.effects.data. In the CreatTraitFiles call, it is given the same name, and passed in as a data table.
	# gen.data is passed in as a string, the name of a file. The file it specifies should have columns of SNP, CLST, FRQ, A1, A2
	gwas.data <- read.table ( gwas.data.filename , header = T , stringsAsFactors = F )
	temp <- unlist ( strsplit ( gen.data , split = "/" ) )
	temp <- temp [ length ( temp ) ]
	
	write ( paste ( gwas.data$SNP , sep = "\n" ) , file = paste ( gwas.data.filename , ".search" , sep = "" ) , ncolumns = 1 )
	my.cmd <- paste ( "perl CodeForRelease/Scripts/sampleSNPs.pl ", gwas.data.filename , ".search " , gen.data , " > ", gsub ( ".txt" , "" , gwas.data.filename) , "." , temp , ".freqs.temp" , sep = "")
	system ( my.cmd )

	# The gen.data file I'm currently using has a number of entries with asterisks after the base ID, which just indicate that the derived allele was not known when Joe Pickrel did the imputation.
	system ( paste ( "sed 's/*//'" , " " , gsub ( ".txt" , "" , gwas.data.filename ) , "." , temp , ".freqs.temp" , " > " , gsub ( ".txt" , "" , gwas.data.filename ) , "." , temp , ".freqs" , sep = "" ) )
	system ( paste ( "rm -f " , gsub ( ".txt" , "" , gwas.data.filename ) , "." , temp , ".freqs.temp" , sep = "" ) )	
	
	
}


RemoveSNPs <- function ( gwas.data.filename , freq.data.filename , output.file ) {
	
	#recover()
	gwas.data <- read.table ( gwas.data.filename , h = T , stringsAsFactors = F )
	freq.data <- read.table ( freq.data.filename , h = T , stringsAsFactors = F )
	have.these <- unique ( freq.data$SNP )
	reduced.gwas.data <- gwas.data [ gwas.data$SNP %in% have.these , ] 
	write.table ( reduced.gwas.data , file = output.file , quote = F , row.names = F , col.names = T )
	
	
	
}

RemovePops <- function ( freq.data.filename , env.var.data.filename ) {
	
	#recover()
	env.var.data <- read.table ( env.var.data.filename , h = T )
	freq.data <- read.table ( freq.data.filename , h = T )
	new.freq.data <- freq.data [ freq.data$CLST %in% env.var.data$CLST , ] 
	write.table ( new.freq.data , file = paste ( freq.data.filename , "." , strsplit ( env.var.data.filename , "/" )[[1]][2] , sep = "" ) , quote = F , row.names = F , col.names = T )
	
	
}