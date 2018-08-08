# DESeq_Custom

Carly D. Kenkel, 
ckenkel@usc.edu

Short Guide 
-----------

1. Put all this into the same directory:
	- script: DESeq.R 
	- table of gene annotations for your sequences: two-column (gene identifier <tab> gene name), 
		tab-delimited, one line per gene. Example given: plob_iso2gene.tab (see annotatingTranscriptomes repository for 		how to generate this file)
	- table of read counts by gene (see tag-based_RNAseq repository for how to generate this file). Example given: 			allcountsPoritesPatchGmapper.txt

2. Make sure you have R installed. The R part requires additional packages, which 
you might need to install prior to running this method. Directions are provided in the .R script

3. Open DESeq.R script; edit the input file names, mark and execute bits of code
separated by blank lines one by one. Follow instructions given as comments in the script.

4. For more information, please refer to the DESeq vignette https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjYvvS4_N3cAhUKrVQKHbPUC0UQFjAAegQIABAC&url=https%3A%2F%2Fbioconductor.org%2Fpackages%2Fdevel%2Fbioc%2Fvignettes%2FDESeq%2Finst%2Fdoc%2FDESeq.pdf&usg=AOvVaw3-roK2pfeegaDc2fuOarq5

5. Don't forget to cite the orignial authors! Within your R console type the following command: 
    citation(DESeq)
   
