

## This README file is made up of 4 components:

	1. Introduction
	2. Overview of folder structure
	3. Overview of data files
	4. Overview of R command scripts 

*Author: Falko Buschke* (`falko.buschke@gmail.com`) *15 September2015*


## 1. Introduction

This archive contains the data and R scripts used for the following study:
	
* Buschke, F.T., Brendonck, L. & Vanschoenwinkel (unpublished). Adding energy gradients and long-distance dispersal to a neutral model improves predictions of Madagascan bird diversity.. 

A written description of the research methodology can be obtained from the manuscript. Any comments or inquiries can be directed to the lead author, Falko Buschke (`falko.buschke@gmail.com`).

While the data and code have been thoroughly tested, we acknowledge that it is by no means perfect. The contents of these files were last verified on 15 September2015.

These files contain pseudo-code in that they are annotated and complete to the extent that they can be used to replicate the results of the aforementioned study. However, they are incomplete in the sense that they cannot recreate the research outputs 100% in their unmodified form. For instance, the simulations were replicated in the study, yet the code does not include this.  


## 2. Overview of folder structure

Once unzipped, the .zip archive - `R_code .zip` - contains one folder (`R_code`) and one sub-folder (`Sensitivity_analysis`)

Unless specified otherwise, all the files are located within the parent directory (`R_code`).

The sub-folder (`Sensitivity_analysis`) contains the data an R scripts for the sensitivity analysis (see Figures 2 in the manuscript).


## 3. Overview of data files

The main directory contains 3 data files in .txt format: `Mada_birds.txt`, `Mada_endemics.txt`, `Mada_grid.txt`.

As well as one sub-folder: `Sensitivity_analysis`


### 3.1. Mada_birds.txt

This dataset contains the community matrix for all birds with complete geographic distribution data in Madagascar, with species as rows and quadrats as columns. There are data for 224 species, which are labelled by their scientific binomial (Species + genus) and ordered alphabetically. The data for 280 quadrats are included (note that labelling begins at 0 and ends at 279) The ID numbers of these quadrates are included in the first column. Each combination of site and species contains a 0 if the species is absent from that quadrat and a 1 if it is present. There are, in total, 62 720 entries of which 37 463 indicate presence (1) and 25 257 indicate absence (0).

### 3.2. Mada_endemics.txt

This dataset contains the community matrix for all birds endemic to Madagascar with complete geographic distribution data, with species as rows and quadrats as columns. There are data for 112 species, which are labelled by their scientific binomial (Species + genus) and ordered alphabetically. The data for 280 quadrats are included (note that labelling begins at 0 and ends at 279) The ID numbers of these quadrates are included in the first column. Each combination of site and species contains a 0 if the species is absent from that quadrat and a 1 if it is present. There are, in total, 31 360 entries of which 16 738 indicate presence (1) and 14 622 indicate absence (0).
	
### 3.3. Mada_grid.txt

This file contains the environmental characteristics of each of the 280 quadrats examined in this study. The dataset is made up of 5 columns and 280 rows (one for each quadrat). The 11 columns represent the following variables:

* 3.1.1. *FID*: This is the Field Identification number for the quadrat. Note that numbering starts at 0 and ends at 279

* 3.1.2. *NDVI*: This represents the Normalised Difference Vegetation Index, an estimate of net primary productivity. These data are from:

	Tucker, C.J., Pinzon, J.E., Brown, M.E., Slayback, D.A., Pak, E.W., Mahoney, R., Vermote, E.F. & Saleous, N.E. (2005) An extended AVHRR 8-km NDVI dataset compatible with MODIS and SPOT Vegetation NDVI Data. *International Journal of Remote Sensing*, **26**, 4485-4498.

* 3.1.3. *Long*: This is the longitude in kilometers (km) 

* 3.1.4. *Latitude*: This is the latitudein kilometers (km)  

* 3.1.5. *Degree*: This is latitude in degrees, for the plot in Figure 1.

	

### 3.4. Subfolder: Sensitivity_analysis

This sub-directory contains 6 data files that contain the output of the sensitivity analysis. 

The following files are in the folder:
* 3.4.1. *Birds_a.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *a* (short-distance dispersal) for all birds in Madagascar.
	
* 3.4.2. *Birds_b.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *b* (long-distance dispersal) for all birds in Madagascar.

* 3.4.3. *Birds_K.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *K* (average habitat capacity) for all birds in Madagascar.

* 3.4.4. *Endemics_a.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *a* (short-distance dispersal) for endemic birds in Madagascar.

* 3.4.5. *Endemics_b.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *b* (long-distance dispersal) for endemic birds in Madagascar.

* 3.4.6. *Endemics_K.txt*: This contains the change in error for the four diversity patterns for incremental changes in parameter *K* (average habitat capacity) for endemic birds in Madagascar.

	

## 4. Overview of R command scripts 

There are 2 fully-annotated files of R command scripts in the `.r` format (these can be opened with any text editor). 

These are: `Neutral_model.r`and `Sensitivity.r`

In order for these command scripts to work, it is important that the data files (Section 3) are all in the same folder (the sensitivity analysis outputs should be in a sub-folder: `Sensitivity analysis`). Then, in each of the scripts, it is essential that the working director be defined explicitly (unless you re using the default directory, in which case you must delete the command line that defines the working directory). 


### 4.1. Neutral_model.r

This contains all the code required to simulate the diversity patterns under neutral assumptions. In order to fully emulate our study, it is necessary to simulate multiple replicates of each type of the model.

**Warning**: The models require a lot of time to run completely!


### 4.2. Sensitivity.r

This contains all the commands necessary to plot the outputs of the sensitivity analyses (Figure 2 in the main manuscript).


