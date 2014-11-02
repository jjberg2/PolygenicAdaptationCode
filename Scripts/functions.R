#setwd("/Users/JeremyBerg/Documents/Academics/CoopLab/Projects/EnvironmentalCorrelations")
#setwd("/Users/JeremyBerg/Documents/Academics/PolygenSel")
PolygenicAdaptationFunction <- function ( gwas.data.file , freqs.file , env.var.data.files , match.pop.file , full.dataset.file , path , match.categories , match.bins  , cov.SNPs.per.cycle = 5000 , cov.cycles = 4 , null.phenos.per.cycle = 1000 , null.cycles = 2 , load.cov.mat = F , PCProj = F , sim.null = T , check.allele.orientation = F ) {
	library ( plyr )
	#recover()
	##### Read in GWAS file, frequencies file, match pop file, and environmental variables files; also sort environmental variable files alphanumerically, create output directory
	dir.create ( path )
	gwas.data <- read.table ( gwas.data.file , header = T , stringsAsFactors = F )
	freqs <- read.table ( freqs.file , header = T , stringsAsFactors = F )
	env.var.data <- lapply ( env.var.data.files , read.table , header = T , stringsAsFactors = F )
	env.var.data <- lapply ( env.var.data , function ( ENV.VAR ) ENV.VAR [ order ( ENV.VAR$CLST ) , ] )
	match.pop <- read.table ( match.pop.file , header = T )
	dir.create ( paste ( path , "/Output" , sep = "" ) )
	
	# sanity checks
	if ( any ( unlist ( lapply ( env.var.data , function ( y ) lapply ( env.var.data , function ( x ) !x$CLST %in% y$CLST ) ) ) ) ) stop ( "Environmental variable datasets differ in the number of populations." )
	if ( any ( !freqs$CLST %in% env.var.data [[ 1 ]]$CLST ) ) stop ( "Frequency dataset has data for populations not present in environmental variable datasets.")
	if ( any ( !env.var.data [[ 1 ]]$CLST %in% freqs$CLST ) ) stop ( "Environmental variable dataset has data for populations not present frequency dataset.")
	if ( any ( !gwas.data$SNP %in% freqs$SNP ) ) stop ( "GWAS dataset contains SNPs that are not in the frequency dataset.")
	if ( any ( !freqs$SNP %in% gwas.data$SNP ) ) stop ( "GWAS dataset contains SNPs that are not in the frequency dataset." )
	
	
	if ( check.allele.orientation == T ) gwas.data <- MatchAlleles ( gwas.data , freqs  )
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
	
	freqs <- freqs [ with ( freqs , order ( SNP , CLST ) ) , ] 
	#recover()
	if ( file.exists ( paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) ) & load.cov.mat == T ) {
		load ( paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) )
	} else {
		uncentered.cov.mat <- SampleCovSNPs ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle = cov.SNPs.per.cycle , cycles = cov.cycles , path = path , full.dataset.file = full.dataset.file , env.var.data = env.var.data )
	}
	T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops - 1 ) ) , times = num.pops - 1 ) , ncol = num.pops , nrow = num.pops - 1 )
	tmp1 <- lapply ( env.var.data , function ( x ) GetOffCovMats ( x , uncentered.cov.mat , gwas.data$EFF , freqs ) )
	regional.off.center.cov.mats <- lapply ( tmp1 , function ( x ) x [[ 1 ]] )
	regional.off.center.T.mats <- lapply ( tmp1 , function ( x ) x [[ 2 ]] )
	tmp2 <- GetOffCovMats ( data.frame ( env.var.data [[ 1 ]] [ , c ( 1 , 2 ) ] , REG = seq_len ( nrow ( env.var.data [[ 1 ]] ) ) ) , uncentered.cov.mat , gwas.data$EFF , freqs )
	individual.off.center.cov.mats <-  tmp2 [[ 1 ]]
	individual.off.center.T.mats <- tmp2 [[ 2 ]]
	F.mat <- T.mat %*% uncentered.cov.mat %*% t ( T.mat )
	F.inv <- solve ( F.mat )
	C.inv <- solve ( t ( chol ( F.mat ) ) )
	#load(file = "CodeForRelease/HeightTest/uncentered.cov.mat.RData")
	epsilon.freqs <- tapply ( freqs$FRQ , freqs$SNP , mean )
	add.gen.var <- sum ( gwas.data$EFF [ order ( gwas.data$SNP ) ]^2  * epsilon.freqs * ( 1 - epsilon.freqs ) )
	
	the.stats <- CalcStats ( 
					freqs = matrix ( freqs$FRQ , nrow = num.pops ) ,
					effects = gwas.data$EFF ,
					env.var.data = env.var.data , 
					var = add.gen.var ,
					uncentered.cov.mat = uncentered.cov.mat ,
					T.mat = T.mat , 
					regional.off.center.cov.mats = regional.off.center.cov.mats , 
					regional.off.center.T.mats = regional.off.center.T.mats , 
					individual.off.center.cov.mats = individual.off.center.cov.mats , 
					individual.off.center.T.mats = individual.off.center.T.mats , 
					path = path , 
					PCProj = PCProj ,
					null = F
				)
	rownames ( the.stats$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( the.stats$reg.Z ) ) )
	colnames ( the.stats$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( the.stats$reg.Z ) ) )
	names ( the.stats$ind.Z ) <- env.var.data[[1]]$CLST
	#recover()
	asymp.p.vals <- list()
	asymp.p.vals$Qx <- 1 - pchisq ( the.stats$Qx , nrow ( T.mat ) )
	asymp.p.vals$reg.Z <- 2 * (1-pnorm ( abs ( the.stats [[ 7 ]] ) ) )
	asymp.p.vals$ind.Z <- 2 * (1-pnorm ( abs ( the.stats$ind.Z ) ) )
	rownames ( asymp.p.vals$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( asymp.p.vals$reg.Z ) ) )
	colnames ( asymp.p.vals$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( asymp.p.vals$reg.Z ) ) )
	names ( asymp.p.vals$ind.Z ) <- env.var.data[[1]]$CLST
	save ( the.stats , file = paste ( path , "/Output/theStats.Robj" , sep = "" ) )
	save ( asymp.p.vals , file = paste ( path , "/Output/asymptotic.pVals.Robj" , sep = "" ) )
	if ( sim.null == F ){
		return ( list ( stats = the.stats , asymptotic.p.vals = asymp.p.vals ) )
	} else {
		p.vals <- list ()
		null.stats <- NullStats ( gwas.data , match.pop , env.var.data , pop.names , uncentered.cov.mat , F.inv , C.inv , T.mat , regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats , bin.names , null.phenos.per.cycle , null.cycles , path , full.dataset.file )
		#recover()
		p.vals$Qx <- sum ( the.stats$Qx < null.stats$Qx ) / length ( null.stats$Qx )
		p.vals$Fst.comp <- sum ( the.stats$Fst.comp < null.stats$Fst.comp ) / length ( null.stats$Fst.comp )
		p.vals$LD.comp <- sum ( the.stats$LD.comp < null.stats$LD.comp ) / length ( null.stats$LD.comp )
		tmp.upper.tail <- rowSums ( null.stats$betas < unlist ( the.stats$betas ) ) / ncol ( null.stats$betas )
		p.vals$betas <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
		p.vals$pearson.rs <- rowSums ( null.stats$pearson.rs^2 > unlist ( the.stats$pearson.rs )^2 ) / ncol ( null.stats$betas )
		tmp.upper.tail <- rowSums ( null.stats$spearman.rhos < unlist ( the.stats$spearman.rhos ) ) / ncol ( null.stats$spearman.rhos )
		p.vals$spearman.rhos <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
		tmp.upper.tail <- mapply ( function ( x , y ) rowSums ( x > y ) / ncol ( x ) , x = null.stats$reg.Z , y = split ( t ( the.stats[[7]] ) , 1:ncol ( the.stats [[ 7 ]]) ) )
		p.vals$reg.Z <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
		rownames ( p.vals$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( p.vals$reg.Z ) ) )
		colnames ( p.vals$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( p.vals$reg.Z ) ) )
		tmp.upper.tail <- rowSums ( null.stats$ind.Z > the.stats$ind.Z ) / ncol ( null.stats$ind.Z )
		p.vals$ind.Z <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
		save ( p.vals , file = paste ( path , "/Output/pVals.Robj" , sep = "" ) )
		save ( null.stats , file = paste ( path , "/Output/nullStats.Robj" , sep = "" ) )
		return ( list ( stats = the.stats , p.vals = p.vals , asymptotic.p.vals = asymp.p.vals ) )
	}
}

GetOffCovMats <- function ( env.var.data , uncentered.cov.mat , effects , freqs ) {
	#recover()
	mat.cols <- list ()
	add.vars.regional <- list ()
	add.vars.individual <- list ()
	k <- 1
	for ( i in sort ( unique ( env.var.data$REG ) ) ) {
		
		num.pops <- nrow ( env.var.data )
		num.center.pops <- sum ( env.var.data$REG != i )
		j <- 1
		mat.cols [[ k ]] <- lapply ( env.var.data$REG , function ( x ) 
										if ( x != i ) { 
											y <- rep ( - 1 / num.center.pops , num.pops ) ; 
											y [ j ] <- ( num.center.pops - 1 ) / num.center.pops ; 
											j <<- j + 1
											return ( y ) 
										} else { 
											y <-  rep ( 0 , num.pops ) ;
											y [ j ] <- 1 ;
											j <<- j + 1 ;
											return ( y ) 
										} 
		)
		epsilons <- with ( subset ( freqs , freqs$CLST %in% env.var.data$CLST [ env.var.data$REG != i ] ) , tapply ( FRQ , SNP , mean ) )
		add.vars.regional [[ k ]] <- 4 * sum ( epsilons * ( 1 - epsilons ) * effects ^2 )
	k <- k + 1
	}
	T.mats <- lapply ( mat.cols , function ( x ) do.call ( cbind , x ) )
	T.mats <- lapply ( T.mats , function ( x ) x [ - sample ( which ( x[1,] != 1 & x[1,] != 0 ) , 1 ) , ] )
	off.center.cov.mats <- lapply ( T.mats , function ( x ) x %*% uncentered.cov.mat %*% t ( x ) )
	return ( list ( off.center.cov.mats , T.mats ) )
}

NullStats <- function ( gwas.data , match.pop , env.var.data , pop.names , uncentered.cov.mat , F.inv , C.inv , T.mat , regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats , bin.names , reps , cycles , path , full.dataset.file ) {
	
	#recover()

	num.pops <- length ( pop.names )
	gwas.cont.table <- table ( gwas.data [ , bin.names ] )
	my.null.bins <- which ( gwas.cont.table != 0 )
	this.many <- reps * gwas.cont.table 
	null.stats <- list ()
	for ( j in 1 : cycles ) {
		sampled.SNPs <- list ()
		effects.list <- list ()
		for ( i in 1 : nrow ( gwas.data ) ) {
			
			# sample a set of <reps> SNPs to match gwas SNP <i> such that they match its properties jointly for all matching categories
			this.snp <- gwas.data [ i , bin.names ]
			in.this.bin <- list ()
			for ( k in 1 : length ( this.snp ) ) {
					in.this.bin [[ k ]] <-  match.pop [ , bin.names [ k ] ] %in% this.snp [[ k ]] 
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
		sampled.null.data <- read.table ( paste ( path , "/null.samples" , j , sep = "" ) , stringsAsFactors = F , h = T  )
		sampled.null.data <- sampled.null.data [ , 1 : 5 ]
		colnames ( sampled.null.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
		sampled.null.data <- sampled.null.data [ sampled.null.data$CLST %in% env.var.data[[1]]$CLST , ]
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
		null.stats [[ j ]] <- mapply ( function ( FREQ , VAR , EFFECTS ) CalcStats ( FREQ , EFFECTS , env.var.data , VAR , uncentered.cov.mat , T.mat , regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats ) , FREQ = null.freqs , VAR = null.vars , EFFECTS = new.effects.list , SIMPLIFY = F )
	}		
	#recover()
	Qx <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 1 ]] ) ) )
	Fst.component <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 2 ]] ) ) )
	LD.component <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 3 ]] ) ) )
	betas <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 4 ]] ) ) ) ) )
	pearson.rs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 5 ]] ) ) ) ) )
	spearman.rhos <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 6 ]] ) ) ) ) )
	reg.Zs <- lapply (  1 : length ( env.var.data ) , function ( u ) do.call ( cbind , lapply ( lapply ( null.stats , function ( w ) lapply ( lapply ( 1 : length ( env.var.data ) , function ( y ) lapply ( w , function ( x ) x [[ 7 ]] [ , y ] ) ) , function ( z ) as.matrix ( do.call ( cbind , z  ) ) ) ) , function ( v ) v [[ u ]] ) ) )
	ind.Zs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) x [[ 8 ]] ) ) ) )
	rownames ( ind.Zs ) <- env.var.data [[ 1 ]]$CLST
	return ( list ( Qx = Qx , Fst.component = Fst.component , LD.component = LD.component , betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Zs = reg.Zs , ind.Zs = ind.Zs ) )

}
PCAProjections <- function ( F.mat , cent.gvs , var , path) {
	#recover()
	my.eig <- eigen ( F.mat * var )
	var.proj <- numeric ( length ( my.eig$values ) )
	for ( i in 1 : length ( my.eig$values ) ) {
		var.proj [ i ] <- ( ( my.eig$vectors [ , i ] %*% cent.gvs ) / sqrt ( my.eig$values [ i ] ) )^2
	}
	dir.create ( paste ( path , "/figs" , sep = "" ) )
	pdf ( paste ( path , "/figs/PCAQx.pdf" , sep = "" ) )
	plot ( var.proj , pch = 20 , xlab = "Principal Component" , ylab = expression ( paste ( "PCA " , Q[X] , sep = "" ) ) )
	dev.off()
	save ( var.proj , file = paste ( path , "/Output/PCVarProj.Rdata" , sep = "" ) )
	save ( my.eig , file = paste ( path , "/Output/PCs.Rdata" , sep = "" ) )
}
CalcStats <- function ( freqs , effects , env.var.data , var , uncentered.cov.mat , T.mat , regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats , path , PCProj = F , null = T ) {
	#recover()
	
	# mean center and drop a pop
	F.mat <- T.mat %*% ( uncentered.cov.mat ) %*% t ( T.mat )
	F.inv <- solve ( F.mat )
	C.inv <- solve ( t ( chol ( F.mat ) ) )
	
	# Qx and components
	contribs <- t ( t ( freqs ) * effects )
	std.contribs <- C.inv %*% T.mat %*% contribs
	adj.cov.mat <- t ( std.contribs ) %*% std.contribs
	Fst.comp <- sum ( diag ( adj.cov.mat ) ) / var
	LD.comp <- sum ( adj.cov.mat [ row ( adj.cov.mat ) != col ( adj.cov.mat ) ] ) / var
	Qx <- Fst.comp + LD.comp 
	
	
	# env var stats
	gvs <- freqs %*% effects
	cent.gvs <- T.mat %*% gvs
	std.gvs <- C.inv %*% cent.gvs
	std.env <- lapply ( env.var.data , function ( x ) C.inv %*% T.mat %*% x$ENV )
	betas <- lapply ( std.env , function ( x ) sum ( std.gvs * x ) / sum ( x^2 ) )
	pearson.rs <- lapply ( std.env , function ( x ) sum ( std.gvs * x ) / sqrt ( sum ( x^2 )*sum ( std.gvs^2 ) ) )
	spearman.rhos <- lapply ( std.env , function ( x ) cor ( std.gvs , x , method = "spearman" ) )
	if ( null == F ) save ( gvs , file = paste ( path , "/Output/genetic.values.Robj" , sep = "" ) )
	if ( PCProj == T ) {
		PCAProjections ( F.mat , cent.gvs , var , path )
	}
	# regional Z scores
	reg.Z.scores <- mapply ( function ( COV , TMAT ) 
		mapply ( function ( THIS.COV , THIS.TMAT ) 
			ZStats ( gvs , var , THIS.COV , THIS.TMAT ) , 
			THIS.COV = COV , THIS.TMAT = TMAT ) ,
		COV = regional.off.center.cov.mats , TMAT = regional.off.center.T.mats
	)
	## individual Z scores
	ind.Z.scores <- mapply ( function ( THIS.COV , THIS.TMAT ) 
		ZStats ( gvs , var , THIS.COV , THIS.TMAT ) , 
		THIS.COV = individual.off.center.cov.mats , THIS.TMAT = individual.off.center.T.mats 
	)
	return ( list ( Qx = Qx , Fst.comp = Fst.comp , LD.component = LD.comp , betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Z = reg.Z.scores , ind.Z = ind.Z.scores ) ) 
}
ZStats <- function ( gvs , var , cov.mat , T.mat ) {
	#recover()
	#cent.gvs <- T.mat %*% gvs
	drop <- which ( apply ( T.mat , 2 , function ( x ) !any ( x <= 1 & x >= 0  ) ) ) # line related to figuring out which population i dropped in earlier steps
	need <- which ( T.mat [ 1 , ] == 1 | T.mat [ 1 , ] == 0 )
	conditioned <- seq_along ( gvs ) [ - c ( need , drop ) ]
	shift.conditioned <- c ( conditioned [ conditioned < drop ] , conditioned [ conditioned > drop ] -1 )
	shift.need <- ifelse ( need > drop , need - 1 , need )
	cond.dist <- condNormal ( gvs [ conditioned ] , rep ( mean ( gvs [ - need ] ) , length ( gvs ) - 1 ) , var * cov.mat , shift.conditioned , shift.need )
	my.Z <- ( sum ( gvs [ need ] ) - sum ( cond.dist$condMean ) ) / sqrt ( sum ( cond.dist$condVar ) )
	return ( my.Z )
	
}
condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
	# Returns conditional mean and variance of x[req.ind] 
	# Given x[given.ind] = x.given
	# where X is multivariate Normal with
	# mean = mu and covariance = sigma
	#recover()
	B <- sigma[req.ind, req.ind]
	C <- sigma[req.ind, given.ind, drop=FALSE]
	D <- sigma[given.ind, given.ind]
	CDinv <- C %*% solve(D)
	cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
	cVar <- B - CDinv %*% t(C)
	return ( list(condMean=cMu, condVar=cVar) )
}

SampleCovSNPs <- function ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle , cycles , path , full.dataset.file , env.var.data ) {
	
	#recover()
	num.pops <- length ( pop.names )
	#T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops ) ) , times = num.pops ) , ncol = num.pops , nrow = num.pops )
	gwas.cont.table <- table ( gwas.data [ , bin.names ] ) 
	my.cov.bins <- which ( gwas.cont.table != 0 )
	total.cov.reps <- ceiling ( 5000 / nrow ( gwas.data ) )
	this.many <- total.cov.reps * gwas.cont.table	
	
	
	
	snp.by.pop <- list ()
	epsilon.cov.snps <- list ()
	var.cov.snps <- list ()
	scaled.snp.by.pop <- list ()
	uncentered.cov.mat <- list ()
	for ( k in 1 : cycles ) {
		sampled.SNPs <- list ()
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
		sampled.cov.data <- read.table ( paste ( path , "/cov.samples" , k , sep = "" ) , stringsAsFactors = F , h = T )
		sampled.cov.data <- sampled.cov.data [ , 1 : 5 ]
		colnames ( sampled.cov.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
		sampled.cov.data <- sampled.cov.data [ sampled.cov.data$CLST %in% env.var.data[[1]]$CLST , ]
		sampled.cov.data$FRQ <- as.numeric ( sampled.cov.data$FRQ )
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
	
	#recover()
	uncentered.cov.mat <- Reduce ( "+" , uncentered.cov.mat)/cycles
	save ( uncentered.cov.mat , file = paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) )
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

