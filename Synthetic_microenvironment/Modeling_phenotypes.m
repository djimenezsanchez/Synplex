function [step] = Modeling_Phenotypes(BasefolderName,LoadfolderName,...
                                   sizeImages,subject_type,subject_number,step,verbose)
% Initialization of variables.
RandIndicesPrev=[];
saveOptimizationLossPheno = [0];
saveOptimizationNumberOfPixelsPheno = [0];

% Load Parameters and cellular neighborhoods mask
load(['Synthetic_microenvironment/parameters/CellPhenotypes_',subject_type,'_',num2str(subject_number),'.mat']);
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat']);
load([BasefolderName,'/Raw/Neighborhoods_GT_',num2str(subject_type),'_nImage_',num2str(subject_number),'_step_',num2str(step),'.mat']);                                             
Contxt_limit = max(PhenoSize)*2;
Phenotype_Abundance = (Phenotype_Abundance'./sum(Phenotype_Abundance'))';

% Initialize mask of cell phenotypes
small = 4;
small_sizeImages = sizeImages./small;
prob = zeros(Ph,1)';prob(end)=1; GT_Pheno=zeros(small_sizeImages(1),small_sizeImages(2));
for Nb_i=1:size(Phenotype_Abundance,1)
    GT_Nb_cut(GT_Nb_cut==Nb_i) = Nb_i;
    small_Neigh = (GT_Nb_cut==Nb_i); small_Neigh = small_Neigh(1:small:end,1:small:end);
    GT_Pheno = GT_Pheno +small_Neigh.*reshape(randsample(Ph,small_sizeImages(1)*small_sizeImages(2),true,Phenotype_Abundance(Nb_i,:)),small_sizeImages);
end
grrd_x = (small_sizeImages(1)-1)/(sizeImages(1)); grrd_y = (small_sizeImages(2)-1)/(sizeImages(2));
[Xq,Yq] = meshgrid((1:grrd_y:small_sizeImages(2)),(1:grrd_x:small_sizeImages(1)));
GT_Pheno = interp2(GT_Pheno,Xq(1:sizeImages(1),1:sizeImages(2)),Yq(1:sizeImages(1),1:sizeImages(2)),'nearest');
Pheno_Cells = zeros(size(GT_Pheno));
Pheno_Nuc = zeros(size(GT_Pheno));

% Calculate stromal orientation surrounding the tumor.
Tumor = GT_Nb_cut==3;
[Pheno_Orientation, Dist_tot] = StromalOrientation(Tumor,max(NeighborhoodMinSize),1,1);
for n_x=10:40:max(NeighborhoodMinSize);
for n_y=10:40:max(NeighborhoodMinSize);
    [Ph_Or, dist] = StromalOrientation(Tumor,max(NeighborhoodMinSize),n_x,n_y);    
    [Ph_Or_4, dist_4] = StromalOrientation(Tumor,max(NeighborhoodMinSize)*8,n_x*4,n_y*4);
    Pheno_Orientation = Pheno_Orientation + Ph_Or + Ph_Or_4/4;
    Dist_tot = Dist_tot + dist + dist_4/4;
end
end
Pheno_Orientation = Pheno_Orientation./(Dist_tot+1e-16);
Pheno_Orientation = Pheno_Orientation-min(Pheno_Orientation(:)); Pheno_Orientation=((Pheno_Orientation/max(Pheno_Orientation(:)))*2)-1;

% Initialize list of painted pixels
excludedIndices = [];
% Contextual Limitation to assign Nb and Phenotypes.
Contxt_limit = max(PhenoSize)*2;
% Obtain indices that are inside the limits of the image.
auxIm = zeros(size(squeeze(GT_Pheno)));
auxIm(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit) = 1; 
Indices =find(auxIm==1);
% Randomice indices to update pixels in a random order
RandIndices = Indices(randperm(length(Indices)));
sizeIM=size(GT_Pheno);

[Imr,Imc,Imz]=size(GT_Pheno);

MaxnNeigh = max([PhenoSize*2]);
excludedIndices = []; %Reinitialize to prevent indices exclueded from COmmunity selection         
cell_number=1; 
 
lastIndex=0; stop_iteration=0;
while ~isempty(RandIndices)

    % Calculate Phenotype as a % of community.
    Perc_Nb_Pheno = zeros(Nb,Ph+1);
    AnatomicalRegion_mean = zeros(Nb,3);
    Pheno_mean = zeros(Ph,3);
    CodificationOfImage=GT_Nb_cut(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit)*100 + GT_Pheno(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit);
    HistoCommPheno = histcounts([CodificationOfImage(:)',[1:Nb*100+Ph]],'BinMethod', 'integers');

    % Measure cell phenotypes in mask
    for iNb=1:Nb
        CodificationOfImage=GT_Nb_cut(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit);
        % Percentage of Anatomical Regions in the entire imagee
        Perc_Nb_Pheno(:,1) = histcounts([CodificationOfImage(:)',[1:Nb]],'BinMethod', 'Integers')';            
        Perc_Nb_Pheno(:,1) = Perc_Nb_Pheno(:,1)./sum(Perc_Nb_Pheno(:,1))*100;
        % Percentage of Phenotypes in each Anatomical Region in the entire imagee
        for iPhen=1:Ph
            Perc_Nb_Pheno(iNb,iPhen+1) = HistoCommPheno(iNb*100+iPhen);
        end
    end
    Perc_Nb_Pheno(:,2:end) = (Perc_Nb_Pheno(:,2:end)./sum(Perc_Nb_Pheno(:,2:end),2))*100;
    Perc_Nb_Pheno(Perc_Nb_Pheno==0)=eps;
    
    % Update interaction-graph of phenotypes in each Nb
    for iNb=1:Nb        
        CellCellInter(:,:,iNb) = CellPhenotypeInteraction(:,:,iNb).*(((Phenotype_Abundance(iNb,:)*100)./Perc_Nb_Pheno(iNb,2:end)).^3);
%         CellCellInter(:,:,iNb) = CellPhenotypeInteraction(:,:,iNb).*(((Phenotype_Abundance(iNb,:)*100)./Perc_Nb_Pheno(iNb,2:end)).^3);
        Diff_Pheno(:,iNb) = mean(abs((Phenotype_Abundance(iNb,:)*100)-Perc_Nb_Pheno(iNb,2:end)));
    end
    
    % Print info
    fprintf(['Loss:',num2str(Diff_Pheno),'_Pixels:',num2str(length(RandIndices)),'\n'])
    saveOptimizationLossPheno = cat(1,saveOptimizationLossPheno,sum(abs(Diff_Pheno)));
    saveOptimizationNumberOfPixelsPheno = cat(1,saveOptimizationNumberOfPixelsPheno,length(RandIndices));
    
    % To lesser number of iterations.
    RandIndices=RandIndices(~ismember(RandIndices,excludedIndices));
    RandIndices(randperm(numel(RandIndices))) = RandIndices;

    RandIndicesPrev = length(RandIndices);
    fprintf([num2str(length(RandIndices)),' '])
    if length(RandIndices)>20               
        Numstep=[20];
    else
        Numstep=length(RandIndices); % Just use the steps necessary to finish the optimization
    end        

    % IF the number of pixels to set werent changed in one step just
    % let the pixels be inserted
    if lastIndex==length(RandIndices) & length(RandIndices)<length(Indices) & stop_iteration>=400
        optimization_is_blocked=true;stop_iteration=380;
    else
        optimization_is_blocked=false; stop_iteration=stop_iteration+1;
    end
    lastIndex=length(RandIndices);

    % Iterate to assign 1.Communities and 2.Phenotypes.
    for i=1:floor(Numstep)
        if ismembc(RandIndices(i),excludedIndices); continue; end            

        % Obtain row/column and check if it is a boundary
        [r,c] = ind2sub([sizeImages(1),sizeImages(2)],RandIndices(i));                                    

        % Calculate cell phenotype based on neighborhood and probability.
        [Neigh_winner,Pheno_winner,Selblock_Orientation,Selblock_Cells] = cell_pheno_neigh(GT_Pheno,r,Contxt_limit,sizeIM,c,...
                                   Pheno_Orientation,Pheno_Cells,GT_Nb_cut,CellCellInter,Ph);                                    

        % Calculate cell morphology
        [Cell,Nucleus] = calc_ellipse(Selblock_Orientation,PhenoSize(Pheno_winner),PhenoEccentricity(Pheno_winner),...
                                      PhenoMorphDeviation(Pheno_winner),RatioNucleousCellSize(Pheno_winner),PhenoPolarity(Pheno_winner));                   

        % Position of cell
        [Ind_Cell,Ind_Cell_Free,R_Cell,C_Cell]=cell_pos(Cell,r,c,GT_Pheno,excludedIndices);
                
        % Position of Nucleus
        [Ind_Nuc] = nuc_pos(Nucleus,Cell,GT_Pheno,r,c);        

        % Are these pixels already been assigned?
        if length(Ind_Cell_Free)>=.85*length(R_Cell) & Pheno_winner<Ph % Paint it if there is more than 85% of available pixels 
            
            % Assign full cell to the tissue 
            GT_Pheno(Ind_Cell_Free) = Pheno_winner;                        
            Pheno_Cells(Ind_Cell_Free) = cell_number;
            Pheno_Nuc(Ind_Nuc) = cell_number; cell_number = cell_number+1;
            excludedIndices = [excludedIndices; Ind_Cell_Free];                 
            
        % Check if there are less than 10% of available pixels.
        elseif (length(Ind_Cell_Free)<=.5*length(R_Cell) | optimization_is_blocked) & length(Ind_Cell_Free)>0
            
        [Cell,Nucleus] = calc_ellipse(Selblock_Orientation,PhenoSize(end),PhenoEccentricity(end),...
                              PhenoMorphDeviation(end),RatioNucleousCellSize(end),PhenoPolarity(end));                           
        [Ind_Cell,Ind_Cell_Free,R_Cell,C_Cell]=cell_pos(Cell,r,c,GT_Pheno,excludedIndices);                        
        [Ind_Nuc] = nuc_pos(Nucleus,Cell,GT_Pheno,r,c);  
            
            % Assign fibers to the available space.
            if Perc_Nb_Pheno(Neigh_winner,end)<Phenotype_Abundance(Neigh_winner,end)*100
                
                GT_Pheno(Ind_Cell_Free) = Ph;                
                Pheno_Cells(Ind_Cell_Free) = 0;
                Pheno_Nuc(Ind_Cell_Free) = 0;                
                excludedIndices = [excludedIndices; r+Imr.*(c-1)]; 
            else
                
                % Assign existing cell to the available space.
                MostCommonCell = GT_Pheno(Ind_Cell); MostCommonCell = mode(MostCommonCell(MostCommonCell~=Ph));
                GT_Pheno(Ind_Cell_Free)= MostCommonCell;
                MostCommonCell = Pheno_Cells(Ind_Cell); MostCommonCell = mode(MostCommonCell(MostCommonCell>0));
                Pheno_Cells(Ind_Cell_Free) = MostCommonCell;                  
                excludedIndices = [excludedIndices; Ind_Cell_Free];                              
            end    
                      
        end
                  
    end               
    % IF true, save tissue neighborhoods at the actual step.
    if verbose==true
        % Save step by step.
        step = step + 1;       
        if mod(step,50)==1      
            Modeling_Phenotypes_Save(BasefolderName,GT_Pheno,Pheno_Nuc,Pheno_Cells,subject_type,subject_number,step,saveOptimizationLossPheno,saveOptimizationNumberOfPixelsPheno)
        end
    end
end
% Save
step = step + 1;
Modeling_Phenotypes_Save(BasefolderName,GT_Pheno,Pheno_Nuc,Pheno_Cells,subject_type,subject_number,step,saveOptimizationLossPheno,saveOptimizationNumberOfPixelsPheno)
             
end


function [Neigh_winner,Pheno_winner,Selblock_Orientation,Selblock_Cells] = cell_pheno_neigh(GT_Pheno,r,Contxt_limit,sizeIM,c,...
                                   Pheno_Orientation,Pheno_Cells,GT_Nb_cut,CellCellInter,Ph)
Selblock = GT_Pheno(max(r-Contxt_limit,Contxt_limit*2):min(r+Contxt_limit,sizeIM(1)-Contxt_limit*2),max(c-Contxt_limit,Contxt_limit*2):min(c+Contxt_limit,sizeIM(2)-Contxt_limit*2));
Selblock_Orientation = Pheno_Orientation(r,c);
Selblock_Cells = Pheno_Cells(max(r-Contxt_limit,Contxt_limit*2):min(r+Contxt_limit,sizeIM(1)-Contxt_limit*2),max(c-Contxt_limit,Contxt_limit*2):min(c+Contxt_limit,sizeIM(2)-Contxt_limit*2));
Selblock_Neigh = GT_Nb_cut(max(r-Contxt_limit,Contxt_limit*2):min(r+Contxt_limit,sizeIM(1)-Contxt_limit*2),max(c-Contxt_limit,Contxt_limit*2):min(c+Contxt_limit,sizeIM(2)-Contxt_limit*2));
CellCellInter_Only_cells = CellCellInter(1:end-1,1:end-1,:);
Selblock_only_Cells = Selblock(Selblock<=(Ph-1));
probComm = mean(CellCellInter_Only_cells(Selblock_only_Cells,:,mode_rand(Selblock_Neigh(:)) )); % Calc Community Mode and Slect Probability.         
[~,Pheno_winner] = max(probComm);
Neigh_winner = mode_rand(Selblock_Neigh(:));
end        


function [Ind_Cell,Ind_Cell_Free,R_Cell,C_Cell]=cell_pos(Cell,r,c,GT_Pheno,excludedIndices)
[R_Cell,C_Cell]=find(Cell==1); R_Cell = R_Cell-size(Cell,1)/2; C_Cell = C_Cell-size(Cell,2)/2;
R_Cell=R_Cell+r; C_Cell=C_Cell+c; %Set the circle to the image.
[Imr,Imc,Imz]=size(GT_Pheno);
Ind_Cell=R_Cell+Imr.*(C_Cell-1)+Imr.*Imc.*(1-1); % Find index position                                                
Ind_Cell_Free = Ind_Cell(~ismember(Ind_Cell,intersect(Ind_Cell,excludedIndices)));
end

function [Ind_Nuc] = nuc_pos(Nucleus,Cell,GT_Pheno,r,c)
[R_Nuc,C_Nuc]=find(Nucleus==1); R_Nuc = R_Nuc-size(Cell,1)/2; C_Nuc = C_Nuc-size(Cell,2)/2;
R_Nuc=R_Nuc+r; C_Nuc=C_Nuc+c; %Set the circle to the image.
[Imr,Imc,Imz]=size(GT_Pheno);
Ind_Nuc=R_Nuc+Imr.*(C_Nuc-1)+Imr.*Imc.*(1-1); % Find index position
Ind_Nuc = Ind_Nuc(Ind_Nuc<numel(GT_Pheno));
end