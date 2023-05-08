
###### DATA_AVAILABILITY ######

    This repository contains 3 folders:
    
    (i)     Phylogenies - contains one script to generate the data and another to build the phylogenies as seen in Figure 2 of the main text.
    (ii)    Results     - contains the statistical data of phylogenies used in the manuscript and the script to generate Figures 3 to 6.
    (iii)   Statistics  - contains a script to generate proper statistical data from a set of phylogenies.

	
    Obs: All R scripts contain the packages required to function correctly.

    




        # PHYLOGENIES 

        Folder contains

		(1) One ".cpp" file:
	
			C++ source code "phylo_patches.cpp"
	
			Requires the GSL scientific library:
		
				Ubuntu systems installation through terminal:
			
					(C++)	~ sudo apt-get install build-essential
					(GSL)	~ sudo apt-get install gsl-bin
					(GSL)	~ sudo apt-get install libgsl-dev
			
				Ubuntu systems compilation through terminal:
		
					~ c++ -o3 phylo_patches.cpp -o [executable_name] -lm -lgsl -lgslcblas
	
			After compilation, we will get an executable.
	
		
			Execution asks terminal input of parameters in the following sequence:
			
				N - number of trees to be generated
				M - mean migration rate
				T - isolation time
			
            For example,
            ./[executable_name] 2 0.08 80
			
            Produces the following "data.txt":
			
				Pop200_Isol80_Mig0.08_MS1_.txt
				Pop200_Isol80_Mig0.08_Temporal1_.txt
				Pop200_Isol80_Mig0.08_MS2_.txt
				Pop200_Isol80_Mig0.08_Temporal2_.txt
				
					
				
		
		(2) One ".R" file:
		
			R script "Trees.R"
			
			Execution asks precisely the same parameters input used for "data.txt" production.
			
            For example:
      
            tree_number = 1,
            isolation_time = 80,
            migration = 0.08,
      
            will generate the complete and extant phylogenies, as in the panels of Figure 2 of the main text.





        # RESULTS 
    
        Folder contains

        (1) 9 files in ".RData" format:

             They account for the mean, confidence interval, and correlation values of the following metrics:
    
            - Number of speciation events
            - Richness
            - Beta-diversity
            - Asymmetry
            - Balance
            - Alpha-value
            - Phylogeny age



        (2) One R script "Figures.R". 
            
            Together with the aforementioned ".RData", it generates the following Figures of the main text:
    
		    FIGURE		DATA
		        3		Comp_(mean,inf and sup).RData
		        4		Ext_(mean,inf and sup).RData
		        5		Over_time80 & Over_time20.RData
		        6		Correlation_values.RData
		        






        # STATISTICS


        Contains the R script "Data_treatment.R". 
        In the file, one will be asked to:
	
        (1) Set the directory to "data.txt" location. 
        (2) Provide the total number of trees and a list of isolation_time and mean_migration_rate parameters.
        
        For example:

        setwd(~/Data/Phylogenies)
        num_trees = 10
        level <- c(80,40,20,0)
        mig   <- c(0,0.002,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04)

        Its execution shall generate the aforementioned 9 statistical data in ".RData" format.
	
