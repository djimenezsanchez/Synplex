function [] = Parameters_Texture(subject_type, subject_number)
% Set parameters for the tissue textue module
load(['Synthetic_microenvironment/parameters/CellPhenotypes_',subject_type,'_',num2str(subject_number),'.mat']);
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat']);

% Marker Names 
Marker_Names = {'DAPI', 'CK', 'T cell'};
Mk = length(Marker_Names);

% Marker expression. Rows: phenotypes. Columns: Markers
Marker_Expression = ones(Ph,Mk).*0;
Marker_Expression(1,1) = .7; % Ph1 expressing Mk1
Marker_Expression(2,1) = .7; % Ph2 expressing Mk1
Marker_Expression(2,2) = .6; % Ph2 expressing Mk2
Marker_Expression(3,1) = .7; % Ph3 expressing Mk1
Marker_Expression(3,3) = .7; % Ph3 expressing Mk3

% nuclear or membrane markers?
Marker_Localization = {'Nuclear','CK','Cytoplasmatic'}; %nuclear or cytoplasmic markers.

% Signal to noise ratio.
SNR = [45,35,45];

% PSF simulation of a microscope. 
NA = 0.45;
Wavelength = [300,461,534]; % In nanometers
gaussFiltImage = [(1/NA)*(Wavelength(1)/1000),(1/NA)*(Wavelength(2)/1000),(1/NA)*(Wavelength(3)/1000)];

% Simulation of leakage between markers. 
gaussFiltLeakage_Right = [0,.2,.2];
gaussFiltLeakage_Left = [.4,0,0];

% PerlinTexture
PerlinPersistence = [0.05,0.03,0.03];

% Initial freq., Final freq.
PerlinFreq = [1,5;  % Mk1
              3,6;  % Mk2
              2,6]; % Mk3

% PerlinTexture_Sparse
PerlinPersistence_Sparse = [.01,0,0]; % Mk1, Mk2, Mk3

% Initial freq., Final freq.
PerlinFreq_Sparse = [1,3;  % Mk1
                     1,3;  % Mk2 
                     1,3];  % Mk3
                     
                 
% Persistence of background noise, Initial freq., and  Final freq.,
BackgroundPerlinNoise = [0.01,2,10];

% Visualize
Parameters_Visualization_Texture(Marker_Expression,Mk,Marker_Names,CellPhenotypeNames)

% Save data. 
save(['Synthetic_microenvironment/parameters/Marker_Expression_',subject_type,'_',num2str(subject_number),'.mat'],...
             'Marker_Localization','Marker_Expression','SNR','gaussFiltImage','gaussFiltLeakage_Right','gaussFiltLeakage_Left',...
             'PerlinPersistence','PerlinPersistence_Sparse','PerlinFreq','PerlinFreq_Sparse','BackgroundPerlinNoise');
end