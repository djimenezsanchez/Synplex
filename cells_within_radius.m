function [vector_interactions] = cells_within_radius(Mask_of_cells, Phenotypes, Cells_Mask,radius)
   
vector_interactions = zeros(1,max(Phenotypes(:)));
Mask_of_cells(isnan(Mask_of_cells))=0;
Cells = regionprops(Mask_of_cells,'Centroid','Area','Perimeter');

[max_x, max_y] = size(Phenotypes);
    
for c = 1:length(Cells)
    if Cells(c).Area>1
        y = Cells(c).Centroid(1);
        x = Cells(c).Centroid(2);
        PH_context = Phenotypes(max(x-radius,1):min(x+radius,max_x),max(y-radius,1):min(y+radius,max_y));    
        Cells_Mask_context = Cells_Mask(max(x-radius,1):min(x+radius,max_x),max(y-radius,1):min(y+radius,max_y));    
        for p = 1:max(Phenotypes(:))
            vector_interactions(p) = vector_interactions(p) + length(unique(Cells_Mask_context.*(PH_context==p)));
        end    
        ph = Phenotypes(round(x),round(y));
        vector_interactions(ph) = vector_interactions(ph) - 1;
    end
end
% vector_interactions = vector_interactions - min(vector_interactions);
vector_interactions = vector_interactions./sum(vector_interactions);

end