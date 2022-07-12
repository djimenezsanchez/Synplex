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
load(['Synthetic_microenvironment/parameters/Marker_Expression_Synthetic_7.mat']);    
MarkerNames = {'Mk1','Mk2','Mk3','Mk4','Mk5','Mk6'};

for subject_type=image_label
    Excel_info = zeros(length(MarkerNames)*(length(CellPhenotypeNames)-2),nImages);
    
for n_image=9
    
    % Load Ground-truth masks
    try    
        Neighborhoods = load([BasefolderName,'/Raw/Neighborhoods_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_1.mat']);
        Cells = load([BasefolderName,'/Raw/Cells_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_2.mat']);
        Phenotypes = load([BasefolderName,'/Raw/Phenotypes_GT_',subject_type{:,1},'_nImage_',num2str(n_image),'_step_2.mat']);
        Markers = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],1);
        Markers(:,:,2) = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],2);
        Markers(:,:,3) = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],3);
        Markers(:,:,4) = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],4);
        Markers(:,:,5) = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],5);
        Markers(:,:,6) = imread([BasefolderName,'/Multiplex/Multiplex_Synthetic_nImage_',num2str(n_image),'.tiff'],6);                        
        Neighborhoods = medfilt2(Neighborhoods.GT_Nb_cut);
        Cells = Cells.Pheno_Cells;
        Phenotypes = Phenotypes.GT_Pheno;
    catch
        continue
    end        
         
    % Quantify cell phenotype abundance 
    Pheno_Marker_Names = {};
    num=1;
    mk_ph=1;
    for Ph=1:length(CellPhenotypeNames)-2
        for Mk=1:length(MarkerNames)
            
            Markers_aux = Markers(:,:,Mk);
               
            % Name of phenotype and marker
            Pheno_Marker_Names{mk_ph} = [CellPhenotypeNames{Ph},' - ',MarkerNames{Mk}];
            
            % Cell phenotype and marker
            Cells_Mask = (Phenotypes==Ph).*Cells;
            Sel_Cells = unique(Cells_Mask );
            Sel_Cells = Sel_Cells(~isnan(Sel_Cells));
            Sel_Cells = Sel_Cells(Sel_Cells>0);
            
            num=1;
                        
            
            for Cl=Sel_Cells'
                
                Excel_info(mk_ph,num) = max(Markers_aux(Cells_Mask==Cl))/2.55;            
                num=num+1;
            end
               
            mk_ph = mk_ph + 1;
        end  
    end    
    
end    
    
    % Save excel gathering information from all
    T = array2table(Excel_info','VariableNames',Pheno_Marker_Names);
    writetable(T,[SavefolderName,'Excel_Cell_Pheno_Marker_',subject_type{:,1},'.xlsx']);
end
    
