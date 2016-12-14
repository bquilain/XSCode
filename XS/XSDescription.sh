#!/bin/zsh +x
#1. Generate the Muon confidence level distributions. If MC/data not created, go into MC/ingmc/Jobs
./GeneratePDFMuCL_Plan -m -i 1 -f 1000
#2. Analyse data and MC w/ these distributions -> Apply XS_CC0pi_Plan (cf MC/ingmc/Jobs and Data for parallal processing)
#3.Apply CC0piSelection. -s select the sample (s=1 for FV selection, s=0 for OV (maybe the latter should be re-tuned because I haven't paid really huge attention to it. If you want real OV selection, use SandSelection_Plan)
./CC0piSelection -s 1
#4. CC0piSelection generate the files needed by Unfolding. Now it can be applied. -d stipulates the data (selected) file to unfold (it can be MC distribution) and -m, the MC to use the efficiency, unfolding matrix etc...
./UnfoldingOptimisation -d files/MCSelected.txt -m files/MCSelected.txt -n 1 -t 20 -s
#n=Number of iteratiibs, s=artificiallaly creates some statistical fluctuations in the -d distribution. The number of toy experiment having these fluctuations is given by -t (20 here)
#Systematic variations can be added in the same way, but instead of -s, one should use -y
#stat variations are done assuming a Poissonian distribution. To control the number of events in the MC and adjust it to the data one, please do it in
