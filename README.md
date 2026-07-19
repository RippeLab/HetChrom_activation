# HetChrom_activation

iMEF CHROMOCENTER PERTURBATION ANALYSIS:
Preprocessing with Fiji:
Images were converted to .tif and maximum projected in Z

Quantification:
Thresholding-based segmentation is based on the functions makeNucMask.R and makeCCMask.R.
4C or 5C R scripts were used depending on whether the source images had 4 or 5 channels (4: e.g. DNA, dCas9, methylation/acetylation; 5: e.g. DNA, dCas9, MSR RNA, AC, p300)
Use scripts in sequence:
- 4C_01_segmentCC uses the two functions to segment nuclei and chromocenters
- 4C_02_curateCC is a semiautomated R script guiding the user through each segmentation so artifacts and erroneous segmentations can be annotated
- 4C_03 + 4C_04 are used for plotting based on the resulting csv files of the previous scripts

- 5C_00_flatfieldcorrection_405: this is a FIJI macro that performs flatfieldcorrection of the DAPI channel (405 nm) to improve segmentation accuracy
- 5C_01_segmentCC: uses the two functions to segment nuclei and chromocenters
- 5C_02_RS-FISH_batch: script using the RS-FISH plugin to call spots. This is used to detect single MSR RNAs.
- 5C_03_curateCC: is a semiautomated R script guiding the user through each segmentation so artifacts and erroneous segmentations can be annotated
- 5C_04_analysis_batch: take the segmentation masks, the spot data and quantifies everything needed using quantNuclei_v01.R
- 5C_05_combinecsv: combines the necessary csvs into one
- 5C_06_plot: used for plotting based on the resulting csv files of the previous scripts

Further annotation can be found within the scripts.


U2OS 2-6-3 REPORTER ARRAY ANALYSIS:
Immunostaining was performed to detect the endogenous proteins RNA Pol II, in the form of its Ser2 (serine 2) and Ser5 (serine 5) modifications. GFP or HaloTag labelling with JaneliaFluor 646 dye indicated the location of the reporter array in the U2OS263 cell lines. Transcription was readout using RNA FISH, using Atto550 labelled oligos.

Preprocessing with Fiji, Cellpose and ilastik:
1. Images were converted to .tif and maximum projected in Z
2. Tif images were split into the individiual channels
3. Cellpose2 [1] was used for segmentation of the nuclei
- 637nm (Halo) and 488nm (GFP) channels were merged and used as input
- Range cell diameter for U2OS 
- Batch processing: ~ 270 diameter  pixels
4. Merged tif images + nuclear masks were imported to ilastik [2], train with 5 images

Quantification:
Individual channel images and Cellpose nuclear masks were used as input in R to quantify image features in regions corresponding to nuclear masks using the custom function quantNuclei_v01.R [3].
RNase controls and background values for each channel was ascertained and then used for thresholding in the co-rec and distal scripts (Reporter_array_analysis_...).

References
Pachitariu, M. & Stringer, C. Cellpose 2.0: how to train your own model. Nat Methods 19, 1634-1641 (2022)
Berg, S., Kutra, D., Kroeger, T. et al. ilastik: interactive machine learning for (bio)image analysis. Nat Methods 16, 1226–1232 (2019)
Frank, L. Perturbing and imaging nuclear compartments to reveal mechanisms of transcription regulation and telomere maintenance. PhD thesis, Heidelberg University (2023)
