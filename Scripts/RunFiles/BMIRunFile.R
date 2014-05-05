source ( "CodeForRelease/Scripts/CreateTraitFile.R")								
source ( "CodeForRelease/Scripts/functions.R")	 
CreateTraitFile ( "Trait_Data/FilesForPaper/BMI.txt" , "Genome_Data/HapMapInHGDP_PositionsAndBValues")
RemoveSNPs ( "Trait_Data/FilesForPaper/BMI.txt" , "Trait_Data/FilesForPaper/BMI.HapMapInHGDP_PositionsAndBValues.freqs" , "Trait_Data/FilesForPaper/newBMI.txt" )	
options ( error = function () {pushover("Error");recover()})
PolygenicAdaptationFunction ( 
									gwas.data.file = "Trait_Data/FilesForPaper/newBMI.txt" , 
									freqs.file = "Trait_Data/FilesForPaper/BMI.HapMapInHGDP_PositionsAndBValues.freqs" , 
									env.var.data.files = list ( "EnvVar/LATS/HGDP_LATS_GLOBAL" ,
																"EnvVar/SUMMERPCS/HGDP_SUMPC1_GLOBAL" ,
																"EnvVar/SUMMERPCS/HGDP_SUMPC2_GLOBAL" ,
																"EnvVar/WINTERPCS/HGDP_WINPC1_GLOBAL" , 
																"EnvVar/WINTERPCS/HGDP_WINPC2_GLOBAL" 
																) , 
									match.pop.file = "Genome_Data/new_French" , 
									full.dataset.file = "Genome_Data/HapMapInHGDP_PositionsAndBValues" , 
									path = "CodeForRelease/BMI" , 
									match.categories = c ( "MAF" , "IMP" , "BVAL" ) ,
									match.bins = list ( seq ( 0 , 0.5 , 0.02 ), c ( 2 ) , seq ( 0 , 1000 , 100 ) ) , 
									cov.SNPs.per.cycle = 5000 , 
									cov.cycles = 4 , 
									null.phenos.per.cycle = 2000 , 
									null.cycles = 5 
									 )
pushover ( "Done with BMI." )