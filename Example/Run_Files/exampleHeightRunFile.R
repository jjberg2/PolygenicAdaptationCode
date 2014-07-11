source ( "~/Documents/Academics/PolygenSel/CodeForRelease/Scripts/CreateTraitFile.R")				
source ( "~/Documents/Academics/PolygenSel/CodeForRelease/Scripts/functions.R")	 
#CreateTraitFile ( "Trait_Data/FilesForPaper/height.txt" , "Genome_Data/HapMapInHGDP_PositionsAndBValues")
#RemoveSNPs ( "Trait_Data/FilesForPaper/height.txt" , "Trait_Data/FilesForPaper/height.HapMapInHGDP_PositionsAndBValues.freqs" , "Trait_Data/FilesForPaper/newheight.txt" )	

options ( error = recover)
PolygenicAdaptationFunction ( 
									gwas.data.file = "Example/Trait_Data/height.txt" , 
									freqs.file = "Example/Trait_Data/height.HapMapInHGDP_PositionsAndBValues.freqs" , 
									env.var.data.files = list ( "Example/EnvVar/HGDP_LATS_GLOBAL" ,
																"Example/EnvVar/HGDP_SUMPC1_GLOBAL" ,
																"Example/EnvVar/HGDP_SUMPC2_GLOBAL" ,
																"Example/EnvVar/HGDP_WINPC1_GLOBAL" , 
																"Example/EnvVar/HGDP_WINPC2_GLOBAL" 
																) , # Note: you can supply as many env.var.data.files concurrently as you want. If only supplying one file it should still be included in a list, e.g. env.var.data.files = list ( "Example/EnvVar/HGDP_LATS_GLOBAL" )
									match.pop.file = "../Genome_Data/French_GLOBAL" , 
									full.dataset.file = "../Genome_Data/HapMapInHGDP_PositionsAndBValues" , 
									path = "Example/Height" , 
									match.categories = c ( "MAF" , "IMP" , "BVAL" ) ,
									match.bins = list ( seq ( 0 , 0.5 , 0.02 ), c ( 2 ) , seq ( 0 , 1000 , 100 ) ) , 
									cov.SNPs.per.cycle = 5000 , 
									cov.cycles = 1 , 
									null.phenos.per.cycle = 1000 , 
									null.cycles = 1 ,
									load.cov.mat = T ,
									sim.null = T
									)