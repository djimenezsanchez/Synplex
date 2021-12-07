function [] = Parameters_Phenotypes(subject_type, subject_number)
% Set parameters for the cell phenotypes module

% Load Information from Neighborhoods.
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat'],'Nb',...
    'Neighborhood_Abundances', 'NeighborhoodMinSize', 'NeighborhoodNames');

CellPhenotypeNames = {'Stromal cells', 'Tumor cells','Stromal-Immune','Lumen','Stromal fibers'};
Ph=length(CellPhenotypeNames); % Number of phenotypes

% Interaction between cell phenotypes. 
CellPhenotypeInteraction = ones(Ph,Ph,Nb).*0.5; %0.5 cells doesnt show interaction.
CellPhenotypeInteraction(2,3,1) = 1; CellPhenotypeInteraction(3,2,1) = 1; % Tumor - stromal-immune attraction
CellPhenotypeInteraction(2,3,2) = 1; CellPhenotypeInteraction(3,2,2) = 1; % Tumor - stromal-immune attraction

% Nucleous radius of cell phenotypes.
PhenoSize = [10,16,13,2,2];

% Ratio cell/nucleous. 
RatioNucleousCellSize = [.6,.9,.7,1,1];

% Eccentricity of cells.
PhenoEccentricity = [.5,.7,.4,1,1];

% Cell complexity. 0: all cells have the same size, 0.2: there is a 20% variation within cell phenotypes.
PhenoMorphDeviation = [.3,.1,.3,0,0];

% Position of nucleous w.r.t. cell.
% if 0 nucleus is centered. if 1 nucleus is on the extremal side.
PhenoPolarity = [0.1,0.2,0.3,0,0];

% Cell Percentage on each neighborhood. rows neighborhoods, columns phenotypes.
Phenotype_Abundance = [60,  0,  10,  0,  30; ... % 'Stroma'
                       0,  95,  0,  0,  5; ... % 'Tumor'
                       0,  0,  0,  100,  0];    % 'Lumen'                                           
Phenotype_Abundance = Phenotype_Abundance./sum(Phenotype_Abundance,2);

% User-defined parameter visualization
Parameters_Visualization_Phenotypes(CellPhenotypeNames,PhenoSize, ...
                                    Ph,Phenotype_Abundance,PhenoEccentricity,Nb,NeighborhoodNames)

% Save Information from Phenotypes
save(['Synthetic_microenvironment/parameters/CellPhenotypes_',subject_type,'_',num2str(subject_number),'.mat'],'CellPhenotypeInteraction','Ph',...
     'CellPhenotypeNames', 'PhenoSize', 'PhenoEccentricity', 'RatioNucleousCellSize', 'Phenotype_Abundance', 'PhenoMorphDeviation','PhenoPolarity');

end


function [] = visualize_real_data_Pheno_Stromal(xyz)
    
    % Histogram of CD8, FoxP3, and CD8-PD1
%     res = histcnd(xyz(:,1),xyz(:,2),xyz(:,4),...
%                     0:max(xyz(:,1))/10:max(xyz(:,1)), ...
%                     0:max(xyz(:,2))/10:max(xyz(:,2)), ...
%                     0:max(xyz(:,4))/10:max(xyz(:,4)));
%     res = histcnd(xyz(:,1),xyz(:,2),xyz(:,4),...
%                     0:33/10:33, ...
%                     0:6.22/10:6.22, ...
%                     0:25/10:25);
    res = histcnd(xyz(:,1),xyz(:,2),xyz(:,4),...
                    0:20/10:20, ...
                    0:20/10:20, ...
                    0:20/10:20);
    res = res./sum(res(:));
    ind = find(res==max(res(:)));    
    [row,col,z] = ind2sub(size(res),ind(1));

    % CD8&FoxP3    
    xyz_h = sum(res,3)';
    x = [1:11]';
    y = [1:11]';    
    [xq,yq] = meshgrid(1:1:11);
    z3 = griddata(x,y,xyz_h,xq,yq,'natural');
    Vq = interp2(z3,4);
    figure; imagesc(Vq(1:length(Vq)/100:length(Vq),1:length(Vq)/100:length(Vq)));
    xlabel('CD8')
    ylabel('FoxP3')
    caxis manual
    caxis([0 0.12]);
    colorbar   
    
   % CD8&CD8-PD0    
    xyz_h = squeeze(sum(res,2));
    x = [1:11]';
    y = [1:11]';    
    [xq,yq] = meshgrid(1:1:11);
    z3 = griddata(x,y,xyz_h,xq,yq,'natural');
    Vq = interp2(z3,4);
    figure; imagesc(Vq(1:length(Vq)/100:length(Vq),1:length(Vq)/100:length(Vq)));
    xlabel('CD8')
    ylabel('CD8-PD1')
    caxis manual
    caxis([0 0.12]);
    colorbar      
    
end

function [] = visualize_real_data_Pheno_Tum(xyz)
    
    % Histogram of CD8, FoxP3, and CD8-FoxP3
%     res = histcnd(xyz(:,1),xyz(:,2),xyz(:,3),...
%                     0:29/10:29, ...
%                     0:1.75/10:1.75, ...
%                     0:2.88/10:2.88);
    res = histcnd(xyz(:,1),xyz(:,2),xyz(:,3),...
                    0:20/10:20, ...
                    0:20/10:20, ...
                    0:20/10:20);
    res = res./sum(res(:));
    ind = find(res==max(res(:)));    
    [row,col,z] = ind2sub(size(res),ind(1));

    % CD8&FoxP3    
    xyz_h = sum(res,3)';
    x = [1:11]';
    y = [1:11]';    
    [xq,yq] = meshgrid(1:1:11);
    z3 = griddata(x,y,xyz_h,xq,yq,'natural');
    Vq = interp2(z3,4);
    figure; imagesc(Vq(1:length(Vq)/100:length(Vq),1:length(Vq)/100:length(Vq)));
    xlabel('CD8')
    ylabel('FoxP3')
    caxis manual
    caxis([0 0.12]);
    colorbar   
    
   % CD8&CD8-FoxP3    
    xyz_h = squeeze(sum(res,2))';
    x = [1:11]';
    y = [1:11]';    
    [xq,yq] = meshgrid(1:1:11);
    z3 = griddata(x,y,xyz_h,xq,yq,'natural');
    Vq = interp2(z3,4);
    figure; imagesc(Vq(1:length(Vq)/100:length(Vq),1:length(Vq)/100:length(Vq)));
    xlabel('CD8')
    ylabel('CD8-FoxP3')
    caxis manual
    caxis([0 0.12]);
    colorbar   
    
    
end
