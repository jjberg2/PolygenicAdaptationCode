source ( "CodeForRelease/Scripts/CreateTraitFile.R")								
source ( "CodeForRelease/Scripts/functions.R")	 
CreateTraitFile ( "Trait_Data/FilesForPaper/T2D.txt" , "Genome_Data/HapMapInHGDP_PositionsAndBValues")
RemoveSNPs ( "Trait_Data/FilesForPaper/T2D.txt" , "Trait_Data/FilesForPaper/T2D.HapMapInHGDP_PositionsAndBValues.freqs" , "Trait_Data/FilesForPaper/newT2D.txt" )	
options ( error = function () {pushover("Error");recover()})
PolygenicAdaptationFunction ( 
									gwas.data.file = "Trait_Data/FilesForPaper/newT2D.txt" , 
									freqs.file = "Trait_Data/FilesForPaper/T2D.HapMapInHGDP_PositionsAndBValues.freqs" , 
									env.var.data.files = list ( "EnvVar/LATS/HGDP_LATS_GLOBAL" ,
																"EnvVar/SUMMERPCS/HGDP_SUMPC1_GLOBAL" ,
																"EnvVar/SUMMERPCS/HGDP_SUMPC2_GLOBAL" ,
																"EnvVar/WINTERPCS/HGDP_WINPC1_GLOBAL" , 
																"EnvVar/WINTERPCS/HGDP_WINPC2_GLOBAL" 
																) , 
									match.pop.file = "Genome_Data/new_French" , 
									full.dataset.file = "Genome_Data/HapMapInHGDP_PositionsAndBValues" , 
									path = "CodeForRelease/T2D" , 
									match.categories = c ( "MAF" , "IMP" , "BVAL" ) ,
									match.bins = list ( seq ( 0 , 0.5 , 0.02 ), c ( 2 ) , seq ( 0 , 1000 , 100 ) ) , 
									cov.SNPs.per.cycle = 5000 , 
									cov.cycles = 4 , 
									null.phenos.per.cycle = 2500 , 
									null.cycles = 4 
									 )
pushover ( "Done with T2D." )