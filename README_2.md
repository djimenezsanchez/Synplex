## Modeling cell phenotypes

To model custom cell phenotypes and their interactions one should edit Parameters_Phenotypes.m located at Synplex/Synthetic_microenvironment/. In this module there are  variables that can be tuned:

### CellPhenotypeNames:
A cell data structure containing names of the cell phenotypes to be simulated e.g.,

`
CellPhenotypeNames = {'Stromal cells', 'Tumor cells','Stromal-Immune','Lumen','Stromal fibers'};
`

- 'Stromal cells', 'Tumor cells', and 'Stromal-Immune' refer to cell phenotype populations 
- 'Lumen' refers to the space without cells that can be found in the backbround of the tissue, or between glands
- 'Stromal fibers' referes to the space between cells. This parameter is useful to simulate different levels of tissue compactness and the extracellular matrix.


### CellPhenotypeInteraction:
A matrix of size "n of phenotypes" x "n of phenotypes" x "n of neighborhoods", specifying the interactions between phenotypes within each of the neighborhoods. 1 is attraction, 0.5 is no interaction, and 0 is repulsion. e.g.,

`
CellPhenotypeInteraction(2,3,1) = 1; CellPhenotypeInteraction(3,2,1) = 1; % Tumor - stromal-immune attraction in the stromal compartment.
`


### PhenoSize
A vector specifying the average size of cell phenotypes. e.g., 

`
PhenoSize = [10,16,13,2,2]; % Nucleus radius in pixels
`

In the case of 'Lumen' and 'Stromal fibers', these values must be defined. 


### RatioNucleousCellSize 
A vector specifying the ratio between the size of the nucleus and the cell. 1 are cells without cytoplasm, 0.5 are cells with equal nucleus and cytoplasm, 0 are cells without nucleus. 

`
RatioNucleousCellSize = [.6,.9,.7,1,1];
`

### PhenoEccentricity 
A vector specifying the eccentricity of cell phenotypes. Ranging from 0 (round cells) to 1 (maximally elongated cells). e.g., 

`
PhenoEccentricity = [.5,.7,.4,1,1];
`

### PhenoMorphDeviation 
A vector specifying the complexity of cell phenotypes. Ranging from 0 (uniform ellipsoidal shape) to 1 (highly irregular shape) for each cell phenotype.

`
PhenoMorphDeviation = [.3,.1,.3,0,0];
`

### PhenoPolarity 
A vector specifying the average cell phenotype polarity. Ranging from 0 (nucleus at the center of the cell) to 1 (nucleus touching the contour of the cell).

`
PhenoPolarity = [0.1,0.2,0.3,0,0];
`

### Phenotype_Abundance
A matrix specifying the abundances of cell phenotypes within each neighborhood. In percentage. e.g., 

`
Phenotype_Abundance = [60,  0,  10,  0,  30; ... % 'Stroma' is composed of 60% of 'stromal cells', 10% 'stromal-immune cells' and 30% 'stromal fibers'
                       0,  95,  0,  0,  5; ... % 'Tumor'
                       0,  0,  0,  100,  0];    % 'Lumen'   
`



