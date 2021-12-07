function [] = Parameters_Texture(subject_type, subject_number)
% Set parameters for the tissue textue module
load(['Synthetic_microenvironment/parameters/CellPhenotypes_',subject_type,'_',num2str(subject_number),'.mat']);
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat']);

% Marker Names 
Marker_Names = {'CK', 'CD8', 'DAPI'};
Mk = length(Marker_Names);

% Marker expression. Rows: phenotypes. Columns: Markers
Marker_Expression = ones(Ph,Mk).*0;
Marker_Expression(1,3) = 1; % Stromal expressing DAPI
Marker_Expression(2,3) = 1; % Tumor expressing DAPI
Marker_Expression(3,3) = 1; % stromal-immune expressing DAPI
Marker_Expression(3,2) = .5; % Stromal-immune expressing CD8
Marker_Expression(5,2) = .05; % Stromal fibers expressing CD8 (sim. non-spcft)
Marker_Expression(2,1) = .8; % Tumor expressing CK

% nuclear or membrane markers?
Marker_Localization = {'CK','Cytoplasmatic','Nuclear'}; %nuclear or cytoplasmic markers.

% Signal to noise ratio.
SNR = [30,30,30];

% PSF simulation of a microscope. 
NA = 0.45;
Wavelength = [461,616,694]; % In nanometers
gaussFiltImage = [(1/NA)*(Wavelength(1)/1000),(1/NA)*(Wavelength(2)/1000),(1/NA)*(Wavelength(3)/1000)];

% Simulation of leakage between markers. 
gaussFiltLeakage = [0.2,.2,0];

% PerlinTexture
PerlinPersistence = [0.2,0.2,0.2]; % CK, CD8, DAPI

% Initial freq., Final freq.
PerlinFreq = [1,7;  % CK
              1,4;  % CD8              
              1,3]; % DAPI 

% PerlinTexture_Sparse
PerlinPersistence_Sparse = [0,0,0.01]; %CK, CD8, DAPI

% Initial freq., Final freq.
PerlinFreq_Sparse = [1,3;  % CK,                     
                     1,3;  % CD8                     
                     1,2]; % DAPI
                 
% Persistence of background noise, Initial freq., and  Final freq.,
BackgroundPerlinNoise = [0.01,2,10];

% Visualize
Parameters_Visualization_Texture(Marker_Expression,Mk,Marker_Names,CellPhenotypeNames)

% Save data. 
save(['Synthetic_microenvironment/parameters/Marker_Expression_',subject_type,'_',num2str(subject_number),'.mat'],...
             'Marker_Localization','Marker_Expression','SNR','gaussFiltImage','gaussFiltLeakage',...
             'PerlinPersistence','PerlinPersistence_Sparse','PerlinFreq','PerlinFreq_Sparse','BackgroundPerlinNoise');
end