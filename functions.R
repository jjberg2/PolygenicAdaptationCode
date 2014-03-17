setwd("/Users/JeremyBerg/Documents/Academics/CoopLab/Projects/EnvironmentalCorrelations")
PolygenicAdaptationFunction <- function ( gwas.data.file , freqs.file , env.var.data.files , match.pop.file , full.dataset.file , path , match.categories , match.bins  , cov.SNPs.per.cycle = 5000 , cov.cycles = 4 , null.phenos.per.cycle = 1000 , null.cycles = 2 ) {
	
	#recover()
	##### Read in GWAS file, frequencies file, match pop file, and environmental variables files
	gwas.data <- read.table ( gwas.data.file , header = T , stringsAsFactors = F )
	freqs <- read.table ( freqs.file , header = T , stringsAsFactors = F )
	env.var.data <- lapply ( env.var.data.files , read.table , header = T , stringsAsFactors = F )
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
	
	
	# get minor allele frequency from effect allele frequency, and flip effects for minor alleles
	gwas.data$MAF <- ifelse ( gwas.data$FRQ > 0.5 , 1 - gwas.data$FRQ , gwas.data$FRQ )	
	gwas.data$MA.EFF <- ifelse ( gwas.data$FRQ > 0.5 , - gwas.data$EFF , gwas.data$EFF )	
	match.pop$MAF <- ifelse ( match.pop$FRQ > 0.5 , 1 - match.pop$FRQ , match.pop$FRQ )	
	
	
	
	# add relevant matching columns that are not in gwas data 
	gwas.data <- cbind ( gwas.data , match.pop [ match ( gwas.data$SNP , match.pop$SNP ) , match.categories [ !match.categories %in% names ( gwas.data ) ] ,drop = FALSE] )

	# assign SNPs to appropriate bins in contingency table
	x <- mapply ( function ( CATS , BINS ) 
				cut ( gwas.data [ , CATS ] , breaks = BINS , include.lowest = T ) ,  
		   CATS = match.categories , BINS = match.bins , SIMPLIFY = FALSE )
	names ( x ) <- bin.names <-  paste ( match.categories , ".BINS" ,sep = "" )
	gwas.data <- cbind ( gwas.data , x )
		
	y <- mapply ( function ( CATS , BINS ) 
				cut ( match.pop [ , CATS ] , breaks = BINS , include.lowest = T ) ,  
		   CATS = match.categories , BINS = match.bins , SIMPLIFY = FALSE )
	names ( y ) <- paste ( match.categories , ".BINS" , sep = "" )
	match.pop <- cbind ( match.pop , y )
	rm ( x , y )
	
	# initialize some handy things
	pop.names <- env.var.data [[ 1 ]]$CLST
	num.pops <- length ( pop.names )
	
	uncentered.cov.mat <- SampleCovSNPs ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle = cov.SNPs.per.cycle , cycles = cov.cycles , path = path , full.dataset.file = full.dataset.file )
	T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops - 1 ) ) , times = num.pops - 1 ) , ncol = num.pops , nrow = num.pops - 1 )
	off.center.cov.mats <- GetOffCovMats ( env.var.data , uncentered.cov.mat )
	F.mat <- T.mat %*% uncentered.cov.mat %*% t ( T.mat )
	F.inv <- solve ( F.mat )
	C.inv <- solve ( t ( chol ( F.mat ) ) )
	#load(file = "CodeForRelease/HeightTest/uncentered.cov.mat.RData")
	epsilon.freqs <- tapply ( freqs$FRQ , freqs$SNP , mean )
	add.gen.var <- 4 * sum ( gwas.data$EFF [ order ( gwas.data$SNP ) ]^2  * epsilon.freqs * ( 1 - epsilon.freqs ) )
	
	NullStats ( gwas.data , match.pop , pop.names , F.inv , C.inv , bin.names , null.phenos.per.cycle , null.cycles , path , full.dataset.file )
	
	
	# CalcCovMat ( cycles = cov.cycles , path = path)
	
}

GetOffCovMats <- function ( env.var.data , uncentered.cov.mat , effects ,  ) {
	recover()
	
	
}

NullStats <- function ( gwas.data , match.pop , pop.names , F.inv , C.inv , bin.names , reps , cycles , path , full.dataset.file ) {
	
	recover()

	num.pops <- length ( pop.names )
	gwas.cont.table <- table ( gwas.data [ , bin.names ] )
	my.null.bins <- which ( gwas.cont.table != 0 )
	this.many <- reps * gwas.cont.table 
	j=1	## for loop over cycles using index j will go in here somewhere
	sampled.SNPs <- list ()
	effects.list <- list ()
	for ( i in 1 : nrow ( gwas.data ) ) {
		
		# sample a set of <reps> SNPs to match gwas SNP <i> such that they match its properties jointly for all matching categories
		this.snp <- gwas.data [ i , bin.names ]
		in.this.bin <- list ()
		for ( k in 1 : length ( this.snp ) ) {
				in.this.bin [[ k ]] <- match.pop [ , bin.names [ k ] ] %in% c(this.snp) [[ k ]] 
		}
		in.this.bin <- do.call ( cbind , in.this.bin )
		matched.SNPs <- match.pop [ rowSums ( in.this.bin ) == length ( bin.names ) , ]$SNP
		sampled.SNPs [[ i ]] <- as.character ( sample ( matched.SNPs , reps , replace = T ) )
		effects.list [[ i ]] <- ifelse ( 
								match.pop [ match ( sampled.SNPs [[ i ]] , match.pop$SNP ) , ]$FRQ >= 0.5, 
								gwas.data$MA.EFF [ i ] ,
								 - gwas.data$MA.EFF [ i ] 
				 )
	}
	all.sampled.SNPs <- unique ( unlist ( sampled.SNPs ) )
	write ( all.sampled.SNPs , file = paste ( path , "/null.SNPs" , j , sep = "" ) , ncolumns = 1 )
	system ( paste ( "Scripts/sampleSNPs.pl " , path , "/null.SNPs" , j , " " , full.dataset.file , " > " , path , "/null.samples" , j , sep = "" ) )
	sampled.null.data <- read.table ( paste ( path , "/null.samples" , j , sep = "" ) , stringsAsFactors = F )
	sampled.null.data <- sampled.null.data [ , 1 : 5 ]
	colnames ( sampled.null.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
	sampled.null.data <- sampled.null.data [ with ( sampled.null.data  , order ( SNP , CLST ) ) , ]
	SNP.table <- lapply ( sampled.SNPs , function ( x ) table ( x ) )
	split.sampled.null.data <- split ( sampled.null.data , sampled.null.data$CLST )
	null.SNP.mat <- do.call ( cbind , sampled.SNPs )
	null.freqs <- lapply ( split.sampled.null.data , function ( this.pop ) 
			apply ( null.SNP.mat , 2 , function ( x )  this.pop [ match ( x , this.pop$SNP ) , "FRQ" ] )
			)
	effects.mat <- do.call ( cbind , effects.list )
	new.effects.list <- alply ( effects.mat ,1 , function ( x ) x )
	null.gvs <- lapply ( null.freqs , function ( x ) rowSums ( x * effects.mat ) )
	mean.null.freqs <- Reduce ( "+" , null.freqs ) / length ( null.freqs )
	null.vars <- rowSums ( effects.mat^2 * mean.null.freqs * ( 1 - mean.null.freqs ) )
	null.freqs <- alply ( aperm ( laply ( null.freqs , function (x ) x ) , c ( 1, 3, 2 ) ) , 3 , function ( y ) y )
	null.stats <- mapply ( function ( FREQ , VAR , EFFECTS ) CalcStats ( FREQ , EFFECTS , env.var.data , VAR , F.inv , C.inv ) , FREQ = null.freqs , VAR = null.vars , EFFECTS = new.effects.list )
	
	
}

CalcStats <- function ( freqs , effects , env.var.data , var , F.inv , C.inv ) {
	recover()
	cent.freqs <- t ( t ( freqs ) - colMeans ( freqs ) ) [ - nrow ( freqs ) , ]
	contribs <- t (  t ( cent.freqs ) * effects )
	std.contribs <- C.inv %*% contribs
	adj.cov.mat <- t ( std.contribs ) %*% std.contribs / var
	Fst.component <- sum ( diag ( adj.cov.mat ) )
	LD.component <- sum ( adj.cov.mat [ row ( adj.cov.mat ) != col ( adj.cov.mat ) ] )
	Qx <- sum ( adj.cov.mat )
	cent.envs <- lapply ( env.var.data , function ( x ) C.inv %*% ( x$ENV - mean ( x$ENV ) ) [ - length ( x$ENV ) ] )
	std.envs <- lapply ( cent.envs , function ( x ) x / sd ( x ) )
	betas <- lapply ( std.envs , function ( x )  lm ( c ( std.contribs ) ~ rep ( x , ncol ( std.contribs ) ) )$coefficients [ 2 ] )
	gvs <- rowSums ( contribs )
	std.gvs <- C.inv %*% gvs 
	pearson.r <- lapply ( std.envs , function ( x ) cor ( x , std.gvs )^2 )
	spearman.rho <- lapply ( std.envs , function ( x ) cor.test ( x , std.gvs , method = c ( "spearman" ) )$estimate )
	
}

SampleCovSNPs <- function ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle , cycles , path , full.dataset.file ) {
	
	#recover()
	num.pops <- length ( pop.names )
	#T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops ) ) , times = num.pops ) , ncol = num.pops , nrow = num.pops )
	gwas.cont.table <- table ( gwas.data [ , bin.names ] ) 
	my.cov.bins <- which ( gwas.cont.table != 0 )
	total.cov.reps <- ceiling ( 5000 / nrow ( gwas.data ) )
	this.many <- total.cov.reps * gwas.cont.table	
	
	
	sampled.SNPs <- list ()
	snp.by.pop <- list ()
	epsilon.cov.snps <- list ()
	var.cov.snps <- list ()
	scaled.snp.by.pop <- list ()
	uncentered.cov.mat <- list ()
	for ( k in 1 : cycles ) {
		j = 1 
		for ( BIN in my.cov.bins ) {
			this.many [ BIN ]
			this.bin <- mapply ( function ( x , y ) x [ y ] ,
			x = dimnames ( gwas.cont.table ) , y = arrayInd ( BIN , dim ( this.many ) ) , SIMPLIFY = FALSE )
			in.this.bin <- list ()
			for ( i in 1 : length ( this.bin ) ) {
				in.this.bin [[ i ]] <- match.pop [ , bin.names [ i ] ] %in% this.bin [ i ] 
			}
			in.this.bin <- do.call ( cbind , in.this.bin )
			matched.SNPs <- match.pop [ rowSums ( in.this.bin ) == length ( bin.names ) , ]$SNP
			sampled.SNPs [[ j ]] <- as.character ( sample ( matched.SNPs , this.many [ BIN ] , replace = T ) )
			j = j + 1 	
		}
		sampled.SNPs.count <- table ( unlist ( sampled.SNPs ) )
		sampled.SNPs <- unlist ( sampled.SNPs )
		write ( unlist ( sampled.SNPs ) , file = paste ( path , "/cov.SNPs" , k , sep = "" ) , ncolumns = 1 )
		system ( paste ( "Scripts/sampleSNPs.pl " , path , "/cov.SNPs" , k , " " , full.dataset.file , " > " , path , "/cov.samples" , k , sep = "" ) )
		sampled.cov.data <- read.table ( paste ( path , "/cov.samples" , k , sep = "" ) , stringsAsFactors = F )
		sampled.cov.data <- sampled.cov.data [ , 1 : 5 ]
		colnames ( sampled.cov.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
		#sampled.SNPs <- read.table ( paste ( path , "/cov.SNPs" , k , sep = "" ) , stringsAsFactors = F ) 
		sampled.cov.data <- sampled.cov.data [ with ( sampled.cov.data , order ( SNP , CLST ) ) , ]
		split.sampled.cov.data <- split ( sampled.cov.data$FRQ , sampled.cov.data$SNP )
		cov.freqs <- mapply ( rep , x = split.sampled.cov.data , times = sampled.SNPs.count )
		snp.by.pop [[ k ]] <- t ( matrix ( unlist ( cov.freqs ) , nrow = num.pops ) )
		epsilon.cov.snps [[ k ]] <- apply ( snp.by.pop [[ k ]] , 1 , mean )
		
		# average of ratios
		var.cov.snps [[ k ]] <- epsilon.cov.snps [[ k ]] * ( 1 - epsilon.cov.snps [[ k ]] )
		scaled.snp.by.pop [[ k ]] <- snp.by.pop [[ k ]] / c ( sqrt ( var.cov.snps [[ k ]] ))
		uncentered.cov.mat [[ k ]] <- cov ( scaled.snp.by.pop [[ k ]] )
	}
	uncentered.cov.mat <- do.call ( "+" , uncentered.cov.mat)/cycles
	return ( uncentered.cov.mat )
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
PolygenicAdaptationFunction ( 
									gwas.data.file = "Trait_Data/europe.height.168" , 
									freqs.file = "Trait_Data/europe.height.HapMapInHGDP_PositionsAndBValues.freqs" , 
									env.var.data.files = "EnvVar/LATS/HGDP_LATS_GLOBAL" , 
									match.pop.file = "Genome_Data/French_GLOBAL" , 
									full.dataset.file = "Genome_Data/HapMapInHGDP_PositionsAndBValues" , 
									path = "CodeForRelease/HeightTest" , 
									match.categories = c ( "MAF" , "IMP" ) ,
									match.bins = list ( seq ( 0 , 0.5 , 0.02 ), c ( 2 ) ) , 
									cov.SNPs.per.cycle = 5000 , 
									cov.cycles = 1 , 
									null.phenos.per.cycle = 1000 , 
									null.cycles = 2 
									 )	
}