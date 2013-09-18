PolygenicAdaptationFunction <- function ( gwas.data.file , freqs.file , env.var.data.files , match.pop.file , full.dataset.file , path ) {
	
	recover()
	##### Read in GWAS file, frequencies file, match pop file, and environmental variables files
	gwas.data <- read.table ( gwas.data.file , header = T , stringsAsFactors = F )
	freqs <- read.table ( freqs.file , header = T , stringsAsFactors = F )
	env.var.data <- lapply ( env.var.data.files , read.table , header = T )
	match.pop <- read.table ( match.pop.file , header = T )
	
	# sanity checks
	if ( any ( unlist ( lapply ( env.var.data , function ( y ) lapply ( env.var.data , function ( x ) !x$CLST %in% y$CLST ) ) ) ) ) stop ( "Environmental variable datasets differ in the number of populations." )
	if ( any ( !freqs$CLST %in% env.var.data [[ 1 ]]$CLST ) ) stop ( "Frequency dataset has data for populations not present in environmental variable datasets.")
	if ( any ( !env.var.data [[ 1 ]]$CLST %in% freqs$CLST ) ) stop ( "Environmental variable dataset has data for populations not present frequency dataset.")
	if ( any ( !gwas.data$SNP %in% freqs$SNP ) ) stop ( "GWAS dataset contains SNPs that are not in the frequency dataset.")
	if ( any ( !freqs$SNP %in% gwas.data$SNP ) ) stop ( "GWAS dataset contains SNPs that are not in the frequency dataset.")
	
	# bit to deal with monomorphic SNPs here, although probably that's a preprocessing step
	
	
	# bit to deal with linked SNPs, although probably also preprocessing
	
	
	# bit to flip alleles so they are the same way round in both the GWAS dataset and freqs dataset	
	# TODO: check on MatchAlleles function and see if it needs to be cleaned up
	gwas.data <- MatchAlleles ( gwas.data , freqs  )
	save ( gwas.data , file = paste ( path , "/GWAS.file.RData" , sep = "" ) )
	
	
	# get minor allele frequency from effect allele frequency, and assign SNPs to appropriate bins in contingency table
	match.pop$MAF <- ifelse ( match.pop$FRQ > 0.5 , 1 - match.pop$FRQ , match.pop$FRQ )	
	
	
	
	
}

MatchAlleles <- function ( gwas.data , assoc.loci.freqs ) {
	
	#protect <- gwas.data
	#protect -> gwas.data
	correct.orientation <- subset ( assoc.loci.freqs , CLST == assoc.loci.freqs$CLST [ 1 ] )
	correct.orientation <- correct.orientation [ order ( correct.orientation$SNP ) , ]
	#gwas.data <- gwas.data [ match ( correct.orientation$SNP , gwas.data$SNP ) , ]
	gwas.data <- gwas.data [ order ( gwas.data$SNP ) , ]
	if ( nrow ( gwas.data ) != nrow ( correct.orientation ) ) {
		stop ( "GWAS dataset and Orientation Matching dataset contain differing numbers of SNPs" )
	}
	
	flip <- gwas.data$A1 == correct.orientation$A2 & gwas.data$A2 == correct.orientation$A1
	dont.flip <- gwas.data$A1 == correct.orientation$A1 & gwas.data$A2 == correct.orientation$A2
	for ( i in 1 : nrow ( gwas.data ) ) {
		if ( flip [ i ] ) {
			gwas.data$A1 [ i ] <- correct.orientation$A1 [ i ]
			gwas.data$A2 [ i ] <- correct.orientation$A2 [ i ]
			gwas.data$EFF [ i ] <- - gwas.data$EFF [ i ]
			gwas.data$FRQ [ i ] <- 1 - gwas.data$FRQ [ i ]
		} else if ( dont.flip [ i ] ) {
			#do nothing
		} else {
			stop ( "Strand Issue")		
		}
	}
	# if ( sum ( flip ) != 0 ) {
		# cat ( "The following allele orientations were flipped:\n" , file = log.file )
		# cat ( c ( gwas.data [ flip , "SNP" ] ) , sep = "\n" , file = log.file )	
	# }
	return ( gwas.data )
}
if ( FALSE ) {
	
PolygenicAdaptationFunction ( 	gwas.data.file = "Trait_Data/europe.height.168" ,
													freqs.file = "Trait_Data/europe.height.HapMapInHGDP_PositionsAndBValues.freqs" ,
													env.var.data.files = list ( "EnvVar/LATS/HGDP_LATS_GLOBAL" , "EnvVar/WINTERPCS/HGDP_WINPC2_GLOBAL" ) , 
													match.pop.file = "Genome_Data/French" , 
													full.dataset.file = "Genome_Data/HapMapInHGDP_PositionsAndBValues" , 
													path = "/Users/JeremyBerg/Documents/Academics/CoopLab/Projects/EnvironmentalCorrelations/CodeForRelease/HeightTest"
													)
	
	
}