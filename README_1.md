## Modeling cellular neighborhoods

To model custom cellular neighborhoods and their interactions one should edit Parameters_Neighborhoods.m located at Synplex/Synthetic_microenvironment/. In this module there are 4 data structures that can be tuned:

### NeighborhoodName:
A cell containing names of the neighborhoods to be simulated e.g., 
`
NeighborhoodNames = {'Stroma','Tumor','Lumen'};
`
