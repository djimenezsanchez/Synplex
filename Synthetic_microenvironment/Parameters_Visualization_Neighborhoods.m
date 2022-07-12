function [] = Parameters_Visualization_Neighborhoods(NeighborhoodNames,Neighborhoods,Neighborhood_Abundances,Nb)
    % Heatmap showing Neighborhood interactions
    figure;
    h = heatmap(NeighborhoodNames,NeighborhoodNames,Neighborhoods,'Colormap',parula,'ColorLimits',[0 1]);
    h.Title = 'Cellular neighborhoods interaction';
    h.XLabel = 'Neighborhoods';
    h.YLabel = 'Neighborhoods'; 

    % Heatmap showing Neighborhood abundances
    figure; 
    b = bar([Neighborhood_Abundances;zeros(1,Nb)]);
    hold on;
    xticklabels({'Neighborhoods'})
    ylabel('Relative neighborhoood abundance');
    colors = double(squeeze(label2rgb(1:Nb,'parula')))/255;
    for i=1:Nb
        b(i).FaceColor = 'flat';
        b(i).CData = colors(i,:);
    end
    legend(NeighborhoodNames)

end