function [] = Parameters_Neighborhoods(subject_type, subject_number)
% Set parameters for the cellular neighborhoods module

% Interactions between neighborhoods.
NeighborhoodNames = {'Stroma','Tumor','Lumen'};
Nb=length(NeighborhoodNames); % Number of cellular neighborhoods.
Neighborhoods_Interaction = ones(Nb,Nb)*0.5; % A value of 0.5 mean that there is no interaction between neighborhoods
Neighborhoods_Interaction(1,2) = 1; Neighborhoods_Interaction(2,1) = 1;  % Stroma and Tumor attraction
Neighborhoods_Interaction(2,3) = 1; Neighborhoods_Interaction(3,2) = 1;  % Tumor and Lumen attraction
Neighborhoods_Interaction(logical(eye(size(Neighborhoods_Interaction,1))))=1; % Helps making consistent neighborhoods.             

% Minimum size of each neighborhood (in pixels)
NeighborhoodMinSize = [40,40,40]; 

% Neighborhood abundance (in %)
Neighborhood_Abundances = [33,33,33]; 
Neighborhood_Abundances = Neighborhood_Abundances/sum(Neighborhood_Abundances);

% Visualize user-defined parameters
Parameters_Visualization_Neighborhoods(NeighborhoodNames,Neighborhoods_Interaction,Neighborhood_Abundances,Nb)

% Save Information from Communities
save(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat'],'Nb','Neighborhoods_Interaction', 'Neighborhood_Abundances', 'NeighborhoodMinSize', 'NeighborhoodNames');

end
