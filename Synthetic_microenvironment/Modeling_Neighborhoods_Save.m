function [] = Modeling_Neighborhoods_Save(BasefolderName,GT_Nb_cut,image_label,n_image,step,saveOptimizationLossNeigh,saveOptimizationNumberOfPixelsNeigh)

save([BasefolderName,'/Raw/Neighborhoods_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.mat'], 'GT_Nb_cut');
imwrite(label2rgb(GT_Nb_cut),[BasefolderName,'/Raw/Neighborhoods_GT_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'])        
figure; plot(saveOptimizationLossNeigh(2:end),'LineWidth',3); ylabel('Tissue neighborhood modeling loss (in %)'); xlabel('Iterations');
print(gcf,[BasefolderName,'/Raw/NeighborhoodLoss_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'],'-dpng','-r300'); 
figure; plot(saveOptimizationNumberOfPixelsNeigh(2:end),'LineWidth',3,'Color','k'); ylabel('Neighborhood unassigned pixels'); xlabel('Iterations');
print(gcf,[BasefolderName,'/Raw/Neighborhoodpixles_',num2str(image_label),'_nImage_',num2str(n_image),'_step_',num2str(step),'.png'],'-dpng','-r300');                 
close all;   

end