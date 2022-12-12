## Modeling cellular neighborhoods

To model custom cellular neighborhoods and their interactions one should edit Parameters_Neighborhoods.m located at Synplex/Synthetic_microenvironment/. In this module there are 4 variables that can be tuned:

### NeighborhoodName:
A cell containing names of the neighborhoods to be simulated e.g., 

`
NeighborhoodNames = {'Stroma','Tumor','Lumen'};
`

### Neighborhoods_Interaction:
A matrix of size "n of neighborhoods" x "n of neighborhoods", specifying the interactions between neighborhoods. 1 is attraction, 0.5 is no interaction, and 0 is repulsion. e.g., 

`
Neighborhoods_Interaction(1,2) = 1; Neighborhoods_Interaction(2,1) = 1;  % Stroma and Tumor attraction
`

### NeighborhoodMinSize:
Neighbors are iteratively assigned to the image mask. To do so, a contextual ROI is used. This value specifies the size of such contextual ROI. A small would result in smaller neighborhoods, and a bigger value in bigger ones.

`
NeighborhoodMinSize = [40,40,40]; % in pixels
`

### Neighborhood_Abundances:
A vector specifying the abundances of each neighbor in the image.

`
Neighborhood_Abundances = [33,33,33] % In percentages
`
