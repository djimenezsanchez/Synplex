function [Orientation, dist_Centr_Tot] = StromalOrientation(Tumor,Contxt_limit,stride_x,stride_y)

step = Contxt_limit;
Orientation = zeros(size(Tumor));
dist_Centr_Tot = zeros(size(Tumor));
% Divide the tumor in crops
for x=stride_x:step:size(Tumor,1)
for y=stride_y:step:size(Tumor,2)    
    
    % Extract different tumor regions
    x_step = x:min(x+step,size(Tumor,1)); y_step = y:min(y+step,size(Tumor,2));
    Tm = bwlabel(Tumor(x_step,y_step)); 
    centroids = regionprops(Tm,'Centroid','EquivDiameter');    
    
    % Analyze each centroid of each tumor
    for c=1:size(centroids,1)
        
        % Calculate pixel distance to centroid
        [X,Y] = meshgrid(1:length(y_step),1:length(x_step));
        X = X-round(centroids(c).Centroid(1));
        Y = Y-round(centroids(c).Centroid(2));              
        dist_Centr = min(centroids(c).EquivDiameter./(sqrt(X.^2+Y.^2)+1e-12),1);
        dist_Centr_Tot(x_step,y_step) = dist_Centr_Tot(x_step,y_step) + dist_Centr;
        Orientation(x_step,y_step) = Orientation(x_step,y_step) + X./(sqrt(X.^2+Y.^2)+1e-12) .* dist_Centr;
        
    end
    
    Orientation(x_step,y_step) = Orientation(x_step,y_step);
        
end
end
Orientation(Orientation<-1) = -1; Orientation(Orientation>1) = 1;
dist_Centr_Tot(dist_Centr_Tot>1) = 1;


