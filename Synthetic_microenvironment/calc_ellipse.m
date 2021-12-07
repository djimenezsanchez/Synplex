function [Cell, Nucleus] = calc_ellipse(Selblock_Orientation,CellSize,CellEccentricity,Morph_dev,Rat_Nuc_Cell,polarity)

% Calculate cell morphology
deg=(Selblock_Orientation+(rand*0.3-0.1))*pi/2;% + rand*0.05;            
x_axis = round(CellSize+rand*CellSize*.5);
y_axis = max(round(CellSize*CellEccentricity+rand*CellSize*.5),1);
size_x = round(CellSize+2+rand*CellSize*.5);
size_y = round(CellSize+2+rand*CellSize*.5);
[Cell,Nucleus] = sparse_ellipse_with_nuc(x_axis,y_axis,deg,size_x,size_y,Morph_dev+(rand*0.1-0.05),Rat_Nuc_Cell+(rand*0.1-0.05),polarity+(rand*0.1-0.05));    
% figure; imshow(Cell+Nucleus,[]);