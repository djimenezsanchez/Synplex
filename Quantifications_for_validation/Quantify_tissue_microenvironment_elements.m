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
    Excel_info = zeros(length(NeighborhoodNames)+length(CellPhenotypeNames),nImages);
    interactivity_Nbgh = zeros(length(NeighborhoodNames),length(NeighborhoodNames),nImages); % Cellular neighborhood interactions
    Interactivity_Matrix_25pixels = zeros(length(CellPhenotypeNames),length(CellPhenotypeNames),nImages); % Cell phenotype interactions
    
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
    
    % Neighborhood abundance    
    for Nb=1:length(NeighborhoodNames)
        Excel_info(Nb,n_image) = (sum(sum(Neighborhoods==Nb))/(sizeImages(1)*sizeImages(2)))*100;
    end
    
    % Cellular neighborhood interactions.
    interactivity_Nbgh = Differences_between_CellNeighb(n_image,Neighborhoods,interactivity_Nbgh);    
    
    % Quantify cell phenotype abundance    
    for Ph=1:length(CellPhenotypeNames)      
        Pos_cells = (Phenotypes==Ph).*Cells;
        % Obtain matrix of neihborhods
        neighbs_this_ph = zeros(size(Neighborhoods));
        c_Nbb = 1;
        for Nbb=(Phenotype_Abundance(:,Ph)>0)'
            neighbs_this_ph = neighbs_this_ph + (Neighborhoods==c_Nbb).*Nbb;
            c_Nbb = c_Nbb + 1;
        end
        
        % Calculate 
        for n_cel=unique(Pos_cells)'
            if jaccard(Pos_cells==n_cel,logical((Pos_cells==n_cel).*neighbs_this_ph))<0.5
                Pos_cells(Pos_cells==n_cel) = 0;
            end
        end
        Excel_info(Nb+Ph,n_image) = length(unique(Pos_cells));
    end  
    
    % Quantify cell phenotypes interactions.
    for Nb=1:length(NeighborhoodNames)        
        for Ph_1=1:max(Phenotypes(:)); 
            Interactivity_Matrix_25pixels(Ph_1,:,n_image) = ...
            cells_within_radius(Cells.*(Neighborhoods==Nb).*(Phenotypes==Ph_1),Phenotypes,Cells,25);
        end        
    end   
    
end
    % Cell neighborhood interactivity
    interactivity_Nbgh = reshape(interactivity_Nbgh,size(interactivity_Nbgh,1)*size(interactivity_Nbgh,2),size(interactivity_Nbgh,3));
    Excel_info = [Excel_info;interactivity_Nbgh];
    Cell_neigh_Names = {};
    num = 1;
    for n_1=1:length(NeighborhoodNames)
    for n_2=1:length(NeighborhoodNames)
        Cell_neigh_Names{num} = [NeighborhoodNames{n_1},' - ',NeighborhoodNames{n_2},' ',num2str(num)];
        num = num + 1;
    end
    end
    
    % Cell phenotype interactivity
    Interactivity_Matrix_25pixels = reshape(Interactivity_Matrix_25pixels,size(Interactivity_Matrix_25pixels,1)*size(Interactivity_Matrix_25pixels,2),size(Interactivity_Matrix_25pixels,3));
    Excel_info = [Excel_info;Interactivity_Matrix_25pixels];
    Cell_Pheno_Names = {};
    num = 1;
    for n_1=1:length(CellPhenotypeNames)
    for n_2=1:length(CellPhenotypeNames)
        Cell_Pheno_Names{num} = [CellPhenotypeNames{n_1},' - ',CellPhenotypeNames{n_2},' ',num2str(num)];
        num = num + 1;
    end
    end
    
    % Save excel gathering information from all
    names = [NeighborhoodNames,CellPhenotypeNames,Cell_neigh_Names, Cell_Pheno_Names];            
    names{18} = 'Background_Lumen';
    T = array2table(Excel_info','VariableNames',names);
    writetable(T,[SavefolderName,'Excel_Cell_Populations_',subject_type{:,1},'.xlsx']);
end
    
