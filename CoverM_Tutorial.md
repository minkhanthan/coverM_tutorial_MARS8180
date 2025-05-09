# CheckM Tutorial 
#### Hello world, today I will be showing you how to use coverM (https://github.com/wwood/CoverM) to assess relative abundance, reads per base and RPKM of MAGs assembled from single nematode. You will use coverM after assembling metagenomes and binning them using an algorithm of your choice. These DNA is sequenced using PacBio HiFi Long read sequencing method. The genomes are assembled using MetaMDBG. The assembled genomes are aligned using MiniMap2. The genomes are binned using MetaBat 


Dear Holly/Alejandro, the script and data is located on the teaching cluster. The full path is as follows:

	/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM


I assume you have permission for it. 


Fastq raw reads:

	/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM/fastq_raw_reads_289 



Binned MAGs from MetaBat binning algorithm: (The MAGs are .fa, CoverM only accepts .fna. I added a loop in the script to convert all .fa into .fna)
	
	/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM/metabat_bin_289 
	
 

Script:

	/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM/scripts

This is the script: 

	#!/bin/sh

	#SBATCH --job-name="coverM"
	#SBATCH --partition=bik_p
	#SBATCH --nodes=1
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task=24
	#SBATCH --mem-per-cpu=2G
	#SBATCH --time=7-00:00:00
	#SBATCH --mail-user=mkh08049@uga.edu
	#SBATCH --mail-type=END,FAIL
	#SBATCH -e coverm.err-%N
	#SBATCH -o coverm.out-%N

	#Path Variables

	module load CoverM/0.4.0
	module load SAMtools
	module load minimap2
	module load FastANI

	RAW=/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM/fastq_raw_reads_289
	
	BINS=/work/mars8180/mkh08049/github_tutorial_desmoscolex_coverM/fastq_raw_reads_289/	m84221_250316_084501_s4.2_DS03Desmoscolex.289_289
	
	FNA_BIN=/home/mkh08049/github_tutorial_desmoscolex_coverM/fna_bin_289
	mkdir -p "$FNA_BIN"
	cp "$BINS"/*.fa "$FNA_BIN"

	for FILE in ${FNA_BIN}/*.fa; do
	SAMPLE=$(basename "$FILE" .fa)
	FNA_FILES="${FNA_BIN}/${SAMPLE}.fna"

	cp "$FILE" "$FNA_FILES"

	coverm genome \
	 --single ${RAW}/m84221_250316_084501_s4.2_DS03Desmoscolex.289_289.fastq \
	 --genome-fasta-files "$FNA_FILES" \
	 --methods relative_abundance reads_per_base rpkm \
	 --min-covered-fraction 0 \
	 --output-format dense \
	 --threads 24 \
	 > "${SAMPLE}_coverm.tsv"
	done

The output should have each MAGs, Relative Abundance, Reads Per Base and RPKM in their respective columns. I have attached the excel sheet with each of those information and more. 
It is in "Desmoscolex_CoverM_289.xlsx" which can be seen in the tab to your side. You will need to convert it into .csv when you upload it onto R. 

Weird thing is, relative abundance is over 100%. It is strange, I assume it is because some of the MAGs have high contamination. I will have to more research on it. 

## We will visualize the Reads Per Base using R

Please install and load the following packages: 

	install.packages("ggplot2")
	install.packages("tidyr")
	install.packages("pheatmap")
	install.packages("RColorBrewer")
	install.packages("dplyr")
	install.packages("scales")
	library(ggplot2)
	library(tidyr) 
	library(pheatmap)
	library(RColorBrewer)
	library(dplyr)
	library(scales)

Realistically you would only need ggplot2, scale, tidyr and dplr but I just added the other in case I wanted to use pheat map but ended up not doing it in this script.

Load the csv

	data <- read.csv("Desmoscolex_CoverM_289.csv")

We only want MAGs and Reads Per Base column from this .csv for our graph

	df <- data.frame(MAG = data$MAG,
                 Reads_Per_Base = data$Reads.Per.Base)

ggplot only likes it long 

	df.long <- df %>% pivot_longer(-MAG) 

plot

	df.long %>% 
  	ggplot(aes(x = MAG, y = name,
            	 label = value, fill = value)) +
 	 geom_tile(color = "grey") + 
 	 scale_fill_continuous(low = "white",
                       	 high = "darkgreen",
                       	 name = "Reads Per Base") + 
 	 scale_y_discrete(limits = rev) +
 	 scale_x_discrete(position = "top") +
 	 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
   	 axis.text.y = element_text(size = 10) ) +
  	labs(x = "MAG", y = NULL)

#### We will plot and visualize Relative Abundance, Reads Per Base and RPKM all together in a heat map. 

Set a new data frame 

	df.2 <- data.frame(MAG = data$MAG,
                  Reads_Per_Base = data$Reads.Per.Base,
                  Relative_Abundance = data$Relative.Abundance,
                  RPKM = data$RPKM)


Make it long 

	df2.long <- df.2 %>% pivot_longer(-MAG)

Plot 

	df2.long %>%
  	group_by(name) %>%
  	mutate(value_scaled = rescale(value)) %>% 
 	ungroup() %>%
  	ggplot(aes(x = MAG, y = name,
             label = value, fill = value_scaled)) +
  	geom_tile(color = "grey") + 
  	scale_fill_continuous(low = "white",
                        high = "darkgreen",
                        name = "Reads Per Base") + 
  	scale_y_discrete(limits = rev) +
  	scale_x_discrete(position = "top") +
  	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10) ) +
 	 labs(x = "MAG", y = NULL)

This should spit out some cool visuals that you could steer the direction of your decisions. 

##Please note, the MAGs are not classified yet taxonomically. My plan is to classify them using GTDB-Tk. 

Thanks for an awesome class. I learned a lot. More than I could imagine. 
















