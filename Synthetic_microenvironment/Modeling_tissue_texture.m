function [] = Modeling_tissue_texture(BasefolderName,LoadfolderName,SavefolderName,subject_type,subject_number,last_step_phenotypes,last_step_neighborhoods,verbose)

load([BasefolderName,'/Raw/Neighborhoods_GT_',num2str(subject_type),'_nImage_',num2str(subject_number),'_step_',num2str(last_step_neighborhoods),'.mat']);                                             
load([BasefolderName,'/Raw/Phenotypes_GT_',num2str(subject_type),'_nImage_',num2str(subject_number),'_step_',num2str(last_step_phenotypes),'.mat']);                                             
load([BasefolderName,'/Raw/Cells_GT_',num2str(subject_type),'_nImage_',num2str(subject_number),'_step_',num2str(last_step_phenotypes),'.mat']);                                             
load([BasefolderName,'/Raw/Nucleus_GT_',num2str(subject_type),'_nImage_',num2str(subject_number),'_step_',num2str(last_step_phenotypes),'.mat']);                                             
load(['Synthetic_microenvironment/parameters/Marker_Expression_',subject_type,'_',num2str(subject_number),'.mat'],...
             'Marker_Localization','Marker_Expression','SNR','gaussFiltImage','gaussFiltLeakage','PerlinPersistence','PerlinPersistence_Sparse','PerlinFreq','PerlinFreq_Sparse','BackgroundPerlinNoise'); 
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat']); 
load(['Synthetic_microenvironment/parameters/CellPhenotypes_',subject_type,'_',num2str(subject_number),'.mat']);
sizeImages = size(GT_Pheno);
M = size(Marker_Expression);
crop_x=1:sizeImages(1);
crop_y=1:sizeImages(2);

% Separate cell masks into a nucleous, cytoplasm, and membrane.
[Nuclear_Mask, Cytoplasm_Mask, Membrane_Mask] = nuc_cyt_mem_cell(sizeImages, Pheno_Cells, GT_Pheno, Pheno_Nuc, M);
figure; imshow(label2rgb(GT_Nb_cut(crop_x,crop_y)),[])
figure;imshow(label2rgb(Pheno_Nuc(crop_x,crop_y)),[]);
figure;imshow(GT_Pheno(crop_x,crop_y),[]);

% Convert phenotypes to fluorescence marker expression.
Im = convert_phenotypes_to_marker_expression(sizeImages,M,Marker_Localization,Marker_Expression,Nuclear_Mask,Cytoplasm_Mask,Membrane_Mask);
to_rgb(Im(crop_x,crop_y,:));

% Add Perlin Noise to the fluorescence image
Im = add_perlinNoise_background(Im,M,sizeImages,BackgroundPerlinNoise,GT_Nb_cut);
to_rgb(Im(crop_x,crop_y,:));
Im = add_Marker_PerlinNoise(Im,M,sizeImages,Marker_Localization,PerlinPersistence,PerlinFreq,PerlinPersistence_Sparse,PerlinFreq_Sparse);
to_rgb(Im(crop_x,crop_y,:));

% Apply psf 
for i=1:size(Im,3);     
    Im(:,:,i) = imgaussfilt(Im(:,:,i),gaussFiltImage(i));     
end
to_rgb(Im(crop_x,crop_y,:));

% image noise and blurring
for mk=1:size(Im,3); Im(:,:,mk) = awgn(Im(:,:,mk),SNR(mk)); end;
to_rgb(Im(crop_x,crop_y,:));
    
% Marker leakage
Im = reshape(Im, [size(Im,1)*size(Im,2),size(Im,3)]);
for i=1:size(Im,2)
    if gaussFiltLeakage(i)>0
        Im(:,i) = imgaussfilt(Im(:,i),gaussFiltLeakage(i));     
    end
end
Im = reshape(Im, [sizeImages(1),sizeImages(2),size(Im,2)]);
to_rgb(Im(crop_x,crop_y,:));

% Save Image multispectral. Image with multispectral values.
filename = [SavefolderName,'Multiplex_',num2str(subject_type),'_nImage_',num2str(subject_number)];  
if exist([filename,'.tiff'], 'file')==2
    delete([filename,'.tiff']);
end
for chan = 1:M(2)
    Im(:,:,chan) = Im(:,:,chan)-min(min(Im(:,:,chan)));
    Im(:,:,chan) = Im(:,:,chan)/max(max(Im(:,:,chan)));
    imwrite(Im(:,:,chan), [filename,'.tiff'], 'writemode', 'append');
end
 
end


function [Nuclear_Mask, Cytoplasm_Mask, Membrane_Mask] = nuc_cyt_mem_cell(sizeImages, Pheno_Cells, GT_Pheno, Pheno_Nuc, M)

Cytoplasm_Mask = zeros(sizeImages);
Membrane_Mask = zeros(sizeImages);
Nuclear_Mask = zeros(sizeImages);
for cell = unique(Pheno_Cells)'
    
   thiscell = Pheno_Cells==cell;
   thiscell_phenotype = mode(GT_Pheno(thiscell));   
   Cytoplasm_Mask(xor(Pheno_Cells==cell,Pheno_Nuc==cell)) = thiscell_phenotype;      
   Membrane_Mask(edge(thiscell)) = thiscell_phenotype;
   Nuclear_Mask(Pheno_Nuc==cell) = thiscell_phenotype;

end
Nuclear_Mask(Nuclear_Mask==0) = M(1)-1; % Zero-values are assigned to background.
Cytoplasm_Mask(Cytoplasm_Mask==0) = M(1)-1; % Zero-values are assigned to background.
Membrane_Mask(Membrane_Mask==0)=M(1)-1; % Zero-values are assigned to background.


end

function [Im] = convert_phenotypes_to_marker_expression(sizeImages,M,Marker_Localization,Marker_Expression,Nuclear_Mask,Cytoplasm_Mask,Membrane_Mask)

Im = zeros(sizeImages(1),sizeImages(2),M(2));
for marker = 1:M(2)
    if strcmp(Marker_Localization{1,marker},'Nuclear')
        Im(:,:,marker) = reshape(Marker_Expression(Nuclear_Mask(:),marker),[sizeImages(1),sizeImages(2)]);
    elseif strcmp(Marker_Localization{1,marker},'Cytoplasmatic')
        Im(:,:,marker) = reshape(Marker_Expression(Cytoplasm_Mask(:),marker),[sizeImages(1),sizeImages(2)]);
    elseif strcmp(Marker_Localization{1,marker},'CK')
        Im(:,:,marker) = reshape(Marker_Expression(Cytoplasm_Mask(:),marker),[sizeImages(1),sizeImages(2)]);        
        Im(:,:,marker) = Im(:,:,marker) + 0.75.*reshape(Marker_Expression(Membrane_Mask(:),marker),[sizeImages(1),sizeImages(2)]);        
    end
end

end

function [Im] = add_perlinNoise_background(Im,M,sizeImages,BackgroundPerlinNoise,GT_Nb_cut)
for mk = 1:M(2)
    GeneralPerlinNoise = Perlin_Noise(sizeImages(1),sizeImages(2),BackgroundPerlinNoise(2):BackgroundPerlinNoise(3),BackgroundPerlinNoise(1));    
    Im(:,:,mk) = Im(:,:,mk)+GeneralPerlinNoise.*0.15.*(GT_Nb_cut<max(GT_Nb_cut(:))); % Add noise to all markers with the exception of the background neighbo
end
end

function [Im] = add_Marker_PerlinNoise(Im,M,sizeImages,Marker_Localization,PerlinPersistence,PerlinFreq,PerlinPersistence_Sparse,PerlinFreq_Sparse)
for mk=1:M(2)
%     if strcmp(Marker_Localization{1,mk},'Cytoplasmatic')
    if PerlinPersistence_Sparse(mk)>0
        Im(:,:,mk) = Im(:,:,mk).*Perlin_Noise_Sparse(sizeImages(1),sizeImages(2),PerlinFreq_Sparse(mk,1):PerlinFreq_Sparse(mk,2),PerlinPersistence_Sparse(mk));    
    end
    Im(:,:,mk) = Im(:,:,mk).*Perlin_Noise(sizeImages(1),sizeImages(2),PerlinFreq(mk,1):PerlinFreq(mk,2),PerlinPersistence(mk));    
%     elseif strcmp(Marker_Localization{1,mk},'CK')
%         Im(:,:,mk) = Im(:,:,mk).*Perlin_Noise(sizeImages(1),sizeImages(2),PerlinFreq(mk,1):PerlinFreq(mk,2),PerlinPersistence(mk));    
%     elseif strcmp(Marker_Localization{1,mk},'Nuclear')
%         Im(:,:,mk) = Im(:,:,mk).*Perlin_Noise(sizeImages(1),sizeImages(2),PerlinFreq(mk,1):PerlinFreq(mk,2),PerlinPersistence(mk));    
%     end    
end
end

