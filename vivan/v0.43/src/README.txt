#===============================================================================
# ViVan - Viral Variance Analysis
#===============================================================================

# Thank you for choosing to use ViVan for your viral population analysis 

# At this time, we recommend you check out our web-server : www.vivanbioinfo.org
# It utilizes our in-house servers which should reduce running times significantly

# If you still wish to run ViVan locally, below are a few instructions on how to go about doing so

# The main script that you will be working with is the completeAnalysis.py script. it is the main script that runs the analysis

#===============================================================================
# PLEASE UPDATE THE PATHES IN THE COMPLETE ANALYSIS SCRIPT BEFORE YOU CONTINUE 
# they are found near the top of the completeAnalysis.py file
# SCRIPTS_DIR should be the main folder where you unpacked the ViVan package
# BWA_PATH should direct to the the bwa folder
#===============================================================================
# Dependancies (these must be installed for the pipeline to work):
# SAMtools is expected to be installed and in the path
# BWA (v0.7.8) is expected to be installed and found inside the ViVan scripts directory
# ea-utils (https://code.google.com/p/ea-utils/) folder must be installed inside the ViVan scripts directory
# Python v>2.6 is expected to be installed with these modules:
# Biopython, scipy and numpy
#===============================================================================

# Usage : python completeAnalysis.py [arguments]
# -f : The configuration file (see example and fill accordingly)
# -c : Add this argument if adapter clipping and quality trimming should are required prior to alignment
# -a : Add this argument if the input is sequence data and requires both alignment and pileup
# -p : Skip alignment and only perform pileup
# -P : Add this if pileup 2 nucleotide rate has already been performed
# -N : Add this argument if sequence reads with unknown (N) alleles should be discarded
# -V : Add this if variant annotation has already been performed
# -M : Use this flag if your only interested in producing variant metrics for all the samples in the configuration file
# -S : Add this if you want only variants that pass strand bias to be analyzed
#===============================================================================

