## Modeling from real data

Twelve sections from formalin-fixed, paraffin-embedded high-grade endometrial carcinomas were stained using a six-color multiplex panel targeting key elements of the immune microenvironment: CD8+ T cells, the transcription factor FoxP3, the bona fide T cell activation marker CD137, the programmed cell death-1 (PD-1), cytokeratin (CK), and the nucleus (DAPI). 

Tissue sections were scanned using a PhenoImager HT, capturing 391 1468×1876×6 image fields. For these tumors, patient-level information was available, the most relevant being the classification of the tumors into two groups: ‘POLE-mutated’ with high levels of immunological activity and good response to treatment (77 image fields from 3 patients), and ‘POLE-WT’ (wild-type) with low levels of immunological activity and poor response to treatment (314 image fields from 9 patients).

### Simulation of cellular neighborhoods from a real multiplex image dataset
The three cellular neighborhoods (tumor, stroma, and background) were simulated, with the background, or lumen, located inside the tumor glands, and the stroma placed surrounding the tumor. 
To simulate this spatial configuration, stroma-tumor, and tumor-lumen neighborhood pairs were set with a high level of interactivity, N_int=1, while the rest of the neighborhood pairs were left with no interactions, N_int=0.5. 
To simulate the abundance of tumor, lumen, and stroma, we used the neighborhood-level QuPath quantifications to create a three-dimensional histogram ρ_3, where each axis specifies the presence of each neighborhood in the image set. 
When the simulator generates a new image, it randomly samples one data point from ρ_3, obtaining three values of abundance for the tumor, stroma, and lumen. 
The data point with the highest probability of being picked (i.e., most common combination of neighborhood abundances in the measured image set) is the one with abundance values of 4%, 22%, and 72% for the lumen, stroma, and tumor, respectively. 

### Simulation of cell phenotypes from a real multiplex image dataset
Fifteen phenotypes were simulated from the QuPath phenotype quantifications. Six of these consisted of “all negative”, CD8+, FoxP3+, CD8+FoxP3+, CD8+PD1+, and CD137+ cells, located in the tumor, and another six displayed the same marker expressions but this time located in the stroma. 
Furthermore, a background phenotype (i.e., absence of cells), and two stromal/fiber phenotypes were defined, one in the tumor and one in the stromal region, representing the extracellular matrix. 
As with the cellular neighborhoods, we created a 15-dimensional histogram (ρ_P15,) from the quantifications of QuPath, to model cell phenotype abundance in the patient cohort. 
To simulate the cell morphology of phenotypes, for each cell phenotype, we estimated the average polarity, eccentricity, complexity, size, as well as the marker localization of the cells seen in real images. 
These values were used to define the configuration of Synplex. We then simulated interactions between our phenotypes based on knowledge acquired from previous studies namely, attraction between CD8+ and FoxP3+ cells, CD8+FoxP3+ and FoxP3+ cells, CD8+ and tumor cells and, finally, CD137+ and tumor cells.

![alt text](https://github.com/djimenezsanchez/Synplex/blob/main/supplfig.jpg)

![alt text](https://github.com/djimenezsanchez/Synplex/blob/main/suppltab.jpg)

### Simulation of marker expression and virtual microscopy from a real multiplex image dataset
To define the image quality parameters of the simulation, we extracted key parameters directly from the acquired images, and microscopy set-up parameters from the optics that were used to generate the images, i.e., PhenoImager HT. The marker expression levels, M_p, were measured on positive areas of two tissues to obtain an experimental average marker expression. To calculate the SNR of each marker, we measured the average expression levels of stained tissue areas (i.e., A_signal) and the background (i.e., A_noise), and these values were then introduced into the following equation SNR_(meas.)=20 log_10⁡(A_signal/A_noise). Moreover, we characterized the spectral leakage between adjacent markers. To this end, we measured the marker expression of positive cells and then calculated each marker’s contribution to spectrally adjacent channels. The leakage was calculated as the marker expression ratio (in %) between target positive cells and their expression values in adjacent spectral channels. 
From the microscopy acquisition setup, we extracted the numerical aperture (NA=0.45) of the objective, and the emission filter spectra wavelengths in nanometers. With this information, we calculated the width of the point spread function using the following equation, M_(PSF-width)=0.61λ/NA, where λ is the peak wavelength of the emission spectra of each marker. In addition to all these parameters, we set the Perlin Noise parameters by estimating the frequency and persistence qualitatively for each marker.

![alt text](https://github.com/djimenezsanchez/Synplex/blob/main/suppltab2.jpg)

