-- new macros 2018/02/08

ShrinkXSFormat:
This macro is used to shrink the XSFormat files to keep only Detected
Events (with at least one reconstruction)
The Flux and POT count histograms are kept unchanged to save the normalization of the file.

GenieReweight:
Reads the neut and genie trees and computes the bin-by-bin Genie/NEUT ratio (true binning)
It also gives the cross section curve for both generators: 2nd output

MergeSciBkgWM:
First, run the WMMC with h2o input, then ch input. This macro merges
the 2 output files selecting (ch in grid) + (h2o not in grid)

CreateRecoErrorFiles:
This macro is used to mimic a +1sigma and -1 sigma global variation of the number of events (normalization error)

DrawReconstructedEventsPlots:
Extension of DrawReconstructedEventsPlots with -W and -r options
(water and ratio cross sections)



-- old macros

CompareCalibrationsPM: 
From MC and Data samples, generate a file format which can be easily compared using CompareDataMC_Draw.c

GeneratePDFMuCL: 
Generate the PDF distribution for the MuCL: ./GeneratePDFMuCL -i intial_file_number (e.g, 1) -f final_file_number (e.g, 100) -m

XS_CC0pi: 
Change the format from INGRID data to XS format.

