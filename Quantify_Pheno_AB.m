% Load ground-truth masks and search for simulated interactions
BasefolderName = 'Images'; check_folder_state(BasefolderName);
LoadfolderName = 'Images/Raw/'; check_folder_state(LoadfolderName);
SavefolderName = 'Images/Multiplex/'; check_folder_state(SavefolderName);
image_label = {'Synthetic'};
nImages = 48;
sizeImages = [1500,1500];
addpath('Quantifications_for_validation')

load(['Synthetic_microenvironment/parameters/CellPhenotypes_Synthetic_1.mat']);    
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_Synthetic_1.mat']);    

for subject_type=image_label
    Excel_info = zeros(length(NeighborhoodNames)*length(CellPhenotypeNames),nImages);
    
for n_image=1:nImages
    
    % Load Ground-truth masks
    try    
        Neighborhoods = load([BasefolderName,'/Raw/Neighborhoods_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_1.mat']);
        Cells = load([BasefolderName,'/Raw/Cells_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_2.mat']);
        Phenotypes = load([BasefolderName,'/Raw/Phenotypes_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_2.mat']);
        Neighborhoods = medfilt2(Neighborhoods.GT_Nb_cut);
        Cells = Cells.Pheno_Cells;
        Phenotypes = Phenotypes.GT_Pheno;
    catch
        continue
    end        
    
     
    % Quantify cell phenotype abundance 
    Pheno_AB_Names = {};
    num=1;
    for Nb=1:length(NeighborhoodNames)
        for Ph=1:length(CellPhenotypeNames)    
            
            % Specify name of cell phenotype ab.
            Pheno_AB_Names{num} = [NeighborhoodNames{Nb},' - ',CellPhenotypeNames{Ph}];
            
            % Obtain individual cells in specific neighborhoood
            Pos_cells = (Phenotypes==Ph).*Cells.*(Neighborhoods==Nb);
            
            Excel_info(num,n_image) = length(unique(Pos_cells));
            
            num=num+1;
                     
        end  
    end    
    
end
    
    
    % Save excel gathering information from all
    T = array2table(Excel_info','VariableNames',Pheno_AB_Names);
    writetable(T,[SavefolderName,'Excel_Cell_Pheno_AB_',subject_type{:,1},'.xlsx']);
end
    
