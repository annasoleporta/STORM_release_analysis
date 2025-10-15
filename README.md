# STORM release analysis

## Citation
This code was used to study protein release from nanocarriers in [1], and it was modified from the original code published in [2]:

[1] A. Solé-Porta, S. Pujals, P. Delcanale, A. Roig. STORM: a Tool to Track Cargo Release from Polymeric Nanocarriers at the Single-Particle Level.

[2] P. Delcanale, L. Albertazzi. DNA-PAINT super-resolution imaging data of surface exposed active sites on particles. *Data in Brief*. **2020**, *30*, 105468. https://doi.org/10.1016/j.dib.2020.105468.

If you use this code in your research or publish any results obtained with it, please cite both papers above.

## Introduction

This code was developed to analyze single-particle data obtained from stochastic optical reconstruction microscopy (STORM). Its main purpose is to detect individual nanocarriers (in our case, PLGA nanocapsules), determine their size and position, and quantify how many cargo molecules (in our case, the protein BSA) remain associated with each one. We used Cyanine 5 and Alexa Fluor 488 for the labeling of the nanocarriers and the protein, respectively. TetraSpeck microspheres were used as fiducial markers for drift correction.

To do this, the code:

1. Reads localization data from three channels: nanocarriers (**REFERENCE**), proteins (**MAIN**), and fiducial markers (**FIDUCIAL**).
2. Identifies the centers of nanocarriers and fiducials using a mean shift clustering algorithm.
3. Filters nanocarriers based on quality criteria such as size, shape, aggregation, and proximity to fiducials.
4. Accurately estimates each particle’s radius.
5. Assigns protein localizations to individual nanocarriers based on distance.

The result is a particle-by-particle quantification that allows the study of nanocarrier features and cargo release.

## ReadCoords3.m
`ReadCoords3.m` reads raw data from SMLM (TXT or CSV) obtained from SMLM analysis in NIS-elements (Nikon N-STORM) or ONI software, and extracts coordinates of interest for further processing. It generates three TXT files of three columns containing XYT coordinates of the **REFERENCE** (named *647* in the script), **MAIN** (named *488* in the script), and **FIDUCIAL** markers channels (named *Fid* in the script) to use for further processing. N-STORM files (TXT) are supposed to be 26-column, and columns 4-5-13 are read as X-Y-T. ONI files (CSV) are supposed to be 11-column, and columns 3-4-2 are read as X-Y-T.

#### Inputs
- Required input parameters:
    - FileName: name of the file with extension, e.g. `STORM_raw_data_example.txt`
    - InputType: denotes the type of file, type `N-STORM` for Nikon software TXT files, or `ONI` for CSV files from ONI.
- Optional input parameters:
    - STORMref: **REFERENCE** channel name for N-STORM data.
    - STORMmain: **MAIN** channel name for N-STORM data.
    - STORMfid: **FIDUCIAL** markers channel name for N-STORM data.
 
#### Outputs
- `XYTref.txt`: X, Y, T(frames) coordinates of localization in the **REFERENCE** channel (in our case, NCs 647), expressed in nm (Ref matrix).
-	`XYTcoordinates.txt`: X, Y, T(frames) coordinates of localization in the **MAIN** channel (in our case, protein 488), expressed in nm (Main matrix).
-	`XYTfid.txt`: X, Y, T(frames) coordinates of localization in the **FIDUCIAL** markers channel (in our case, TetraSpeck 561), expressed in nm (Fid matrix).

The command line to run the code would be:

    ReadCoords3('STORM_raw_data_example.txt', 'N-STORM', STORMref='405/647', STORMmain='405/488', STORMfid='Bead Drift Correction')


## Cluster2_quantification.m
`Cluster2_quantification.m` performs mean-shift clustering of XYT coordinates of **REFERENCE** and **FIDUCIAL** channels to roughly identify centers of valid particles and fiducial markers. Then, the localizations in **MAIN** channel within a defined distance from centers are stored. 

#### Inputs
- Required input parameters:
    - FileNameMain: name of the file with XYT coords of **MAIN** channel with extension, e.g. `XYTcoordinates.txt`
    - FileNameRef: name of the file with XYT coords of **REFERENCE** channel with extension, e.g. `XYTref.txt`
    - FileNameFid: name of the file with XYT coords of **FIDUCIAL** channel with extension, e.g. `XYTfid.txt`
    - Bandwidth: parameter for clustering, in nm (half radius of particle).
    - MinPts: minimum number of localizations in a cluster.
    - MaxDiam: maximum diameter of clusters in reference channel (longest axis), in nm.
    - FactorMaxDist: factor that multiplies the radius of each NC and defines the maximum distance between cluster center in reference channel and localization in main channel to be considered attached.
 - Optional input parameters:
    - Elong: maximum ellipse elongation allowed, default=10 (long axis:short axis ratio).
    - ScaleFactor: scale factor in ellipse fit, default=1.0 (multiplicative scale factor used in ellipse fitting).
    - MinClustDist: minimum distance between clusters to be considered isolated, in nm, default=300 nm (if distance > MinClustDist, clusters are isolated).
    - DistFidRef: threshold distance to label a reference cluster as a fiducial marker, in nm, default=100 nm (if distance < DistFidRef, the identified cluster in the reference channel is actually a fiducial marker).

### Outputs
- *Fig1*: Plot of the localizations with the results of cluster filtering (clusters having a number of localizations above the threshold MinPts are fitted with an ellipse model). Magenta: selected clusters; Others: discarded (cyan: elongated; yellow: aggregated; blue: fiducial marker).
- *Fig2*: Plot of the valid particles with associated main localizations (protein). Red circle indicates the nanocarrier size, while green circle contains the protein localizations associated with the particle.
- *Fig3*: Histogram showing the number of main (protein) localizations per particle.

Imagine that we choose the following parameters: Bandwidth=50, MinPts=20, Maxdiameter=200, FactorMaxDist=1.3, Elong=10 (default), ScaleFactor=1 (default), MinClustDist=50, DistFidRef=100. The command line to run the could would be:

    Cluster2_quantification('XYTcoordinates.txt', 'XYTref.txt', 'XYTfid.txt', 50, 20, 200, 1.3, MinClustDist=50, DistFidRef=100)
