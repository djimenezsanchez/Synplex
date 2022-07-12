function [] = Modeling_Phenotypes_Save(BasefolderName,GT_Pheno,Pheno_Nuc,Pheno_Cells,image_label,n_image,step,saveOptimizationLossPheno,saveOptimizationNumberOfPixelsPheno)

% Save Cells
save([BasefolderName,'/Raw/Cells_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.mat'], 'Pheno_Cells');
imwrite(Pheno_Cells,[BasefolderName,'/Raw/Cells_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'])        

% Save Nucleus
save([BasefolderName,'/Raw/Nucleus_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.mat'], 'Pheno_Nuc');
imwrite(Pheno_Nuc,[BasefolderName,'/Raw/Nucleus_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'])        

% Save phenotypes
save([BasefolderName,'/Raw/Phenotypes_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.mat'], 'GT_Pheno');
imwrite(GT_Pheno,[BasefolderName,'/Raw/Phenotypes_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'])        

figure; plot(saveOptimizationLossPheno(2:end),'LineWidth',3); ylabel('Cell phenotype modeling loss (in %)'); xlabel('Iterations');
print(gcf,[BasefolderName,'/Raw/PhenotypeLoss_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'],'-dpng','-r300'); 

figure; plot(saveOptimizationNumberOfPixelsPheno(2:end),'LineWidth',3,'Color','k'); ylabel('Phenotype unassigned pixels'); xlabel('Iterations');
print(gcf,[BasefolderName,'/Raw/Phenotypepixles_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'],'-dpng','-r300');                 
close all;   

end