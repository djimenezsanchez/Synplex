clear all;
close all;

%% Folder where synthetic images are stored
BasefolderName = 'Images'; check_folder_state(BasefolderName);  
LoadfolderName = 'Images/Raw/'; check_folder_state(LoadfolderName); % Folder of Ground-truth masks
SavefolderName = 'Images/Multiplex/'; check_folder_state(SavefolderName); % Folder of Multiplex images

%% Parameters of Synplex
nImages = 1;    % Number of images that will be generated
sizeImages = [500,500]; % x and y dimensions of the generated images
image_label={'Synthetic'}; % Types of images
addpath('Synthetic_microenvironment')
verbose = 0; % Show optimization info.

%% Synplex
for subject_type=image_label
for n_image=1:nImages
    
    % Synplex user-defined parameters
    Parameters_Neighborhoods(subject_type{1},n_image) % Creates Synthetic_microenvironment/parameters/CellularNeighborhoods.mat
    Parameters_Phenotypes(subject_type{1},n_image) % Creates Synthetic_microenvironment/parameters/CellPhenotypes.mat        
    Parameters_Texture(subject_type{1},n_image) % Creates Synthetic_microenvironment/parameters/TissueTexture.mat
 
    % Modeling neighborhoods
    last_step_neighborhods = Modeling_Neighborhoods(BasefolderName,LoadfolderName,...
                                sizeImages,subject_type{1},n_image,verbose);  

    % Modeling cell phenotypes
    last_step_phenotypes = Modeling_phenotypes(BasefolderName,LoadfolderName,...
                             sizeImages,subject_type{1},n_image,last_step_neighborhods,verbose);

    % Tisuue texture and virtual microscope
    Modeling_tissue_texture(BasefolderName,LoadfolderName,SavefolderName,subject_type{1},...
                            n_image,last_step_phenotypes,last_step_neighborhods,verbose);                       
end
end
