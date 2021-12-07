function [] = Parameters_Visualization_Texture(Marker_Expression,Mk,Marker_Names,CellPhenotypeNames)

    figure; 
    b = bar(categorical(CellPhenotypeNames),Marker_Expression);    
    hold on;
    ylabel('Marker Expression');
    colors = [0 0 1;
              0 1 0;
              1 0 0;
              1 1 0;
              1 0 1;
              0 1 1;
              0 0 0];
    for i=1:Mk    
        b(i).FaceColor = 'flat';
        b(i).CData = colors(i,:);
    end
    legend(Marker_Names)

end