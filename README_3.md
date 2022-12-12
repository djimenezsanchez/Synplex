## Cell expression and virtual microscopy
To model marker expression intensity and simulate the effect of the acquisition by a fluorescence microscope one should edit Parameters_Texture.m located at Synplex/Synthetic_microenvironment/. In this module there are variables that can be tuned:

### Marker_Names
A cell data structure containing names of the markers to be simulated e.g.,

`
Marker_Names = {'DAPI', 'CK', 'T cell'};
`

- 'DAPI' stains the cell nucleus
- 'CK' stains epithelial cells (tumor cells)
- 'T cells' stains t-cell cytoplasm.


### Marker_Expression

A matrix specifying the average marker expression 'n of phenotypes' x 'n of markers', where each row specifies the intensity of the expected markers for one cell phenotype. 

`
Marker_Expression(1,1) = .7; % Ph1 expressing Mk1
`


### Marker_Localization 
A cell structure specifying the localization of markers in the cell. 'Nuclear' if it is expressed in the nuclei, 'Cytoplasmatic' if it is expressed in the cytoplasm, and 'CK' if it is expressed in both compartments. e.g., 

`
Marker_Localization = {'Nuclear','CK','Cytoplasmatic'}; %nuclear or cytoplasmic markers.
`


### SNR
A vector specifying the average Signal-to-Noise Ratio (SNR) for each marker. In dB. e.g., 

`
SNR = [45,35,45];
`


### NA and Wavelength
These values together calculate point-spread-function (PSF) of the simulated microscope. NA is a value specifying the numerical aperture of the system. Wavelenght is a vector specifying wavelenght at which the signal is captured, in nanometers.


### gaussFiltLeakage_Right and gaussFiltLeakage_Left
Two vectors specifying the marker expression leakage between consecutive markers. In gaussFiltLeakage_Right each marker is its own expression and a percentage of value of its left. This is why the most-left value have a value equals zero. 

`
gaussFiltLeakage_Right = [0,.2,.2];
gaussFiltLeakage_Left = [.4,0,0];
`

### PerlinPersistence and PerlinFreq 
Vectors specifying the texture of each of the markers. Persistence is the perceived intensity of the texture. PerlinFreq is the frequencies used to generate the noise.

`
PerlinFreq = [1,5;  % Mk1
              3,6;  % Mk2
              2,6]; % Mk3
`


