function [Ph2_in_Nb5_RelAb,Ph2_in_Nb5_Cells,Ph4_in_Nb5_RelAb,Ph4_in_Nb5_Cells, ...
    Ph7_in_Nb5_RelAb,Ph7_in_Nb5_Cells] = Differences_phenotype_Ab(Phenotypes, ...
    Phenotypes_eroded,Neighborhoods,iImages,Ph2_in_Nb5_RelAb,Ph2_in_Nb5_Cells, ...
    Ph4_in_Nb5_RelAb,Ph4_in_Nb5_Cells,Ph7_in_Nb5_RelAb,Ph7_in_Nb5_Cells)

% Background 
Background_in_Nb5 = (Neighborhoods==5).*(Phenotypes==9);
% Ph4_Nb5
Ph4_in_Nb5 = (Neighborhoods==5).*(Phenotypes==4);    
Ph4_in_Nb5_RelAb(iImages) = (sum(Ph4_in_Nb5(:))/(sum(sum(Neighborhoods==5))-sum(sum(Background_in_Nb5))))*100;
Ph4_in_Nb5_Cells(iImages) = length(regionprops(bwlabel((Neighborhoods==5).*(Phenotypes_eroded==4)),'Centroid'));
% Ph2_Nb5
Ph2_in_Nb5 = (Neighborhoods==5).*(Phenotypes==2);    
Ph2_in_Nb5_RelAb(iImages) = (sum(Ph2_in_Nb5(:))/(sum(sum(Neighborhoods==5))-sum(sum(Background_in_Nb5))))*100;
Ph2_in_Nb5_Cells(iImages) = length(regionprops(bwlabel((Neighborhoods==5).*(Phenotypes_eroded==2)),'Centroid'));
% Ph7_Nb5
Ph7_in_Nb5 = (Neighborhoods==5).*(Phenotypes==7);    
Ph7_in_Nb5_RelAb(iImages) = (sum(Ph7_in_Nb5(:))/(sum(sum(Neighborhoods==5))-sum(sum(Background_in_Nb5))))*100;
Ph7_in_Nb5_Cells(iImages) = length(regionprops(bwlabel((Neighborhoods==5).*(Phenotypes_eroded==7)),'Centroid'));

% Boxplot comparing
[h,p4_2] = ttest(Ph2_in_Nb5_Cells,Ph4_in_Nb5_Cells);
[h,p4_7] = ttest(Ph7_in_Nb5_Cells,Ph4_in_Nb5_Cells);
[h,p2_7] = ttest(Ph7_in_Nb5_Cells,Ph2_in_Nb5_Cells);
figure; boxplot([Ph2_in_Nb5_Cells,Ph4_in_Nb5_Cells,Ph7_in_Nb5_Cells],'Labels',{'Ph2','Ph4','Ph7'}); title(strcat('p4_2=',string(p4_2),'p4_7=',string(p4_7),'p2_7=',string(p2_7)));
figure; boxplot([Ph2_in_Nb5_RelAb,Ph4_in_Nb5_RelAb,Ph7_in_Nb5_RelAb],'Labels',{'Ph2','Ph4','Ph7'}); title(strcat('p4_2=',string(p4_2),'p4_7=',string(p4_7),'p2_7=',string(p2_7)));

end