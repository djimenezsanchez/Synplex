function [] = Parameters_Visualization_Phenotypes(CellPhenotypeNames,PhenoSize, CellPhenotypeInteraction,F,...
                                                Phenotype_Abundance,PhenoEccentricity, PhenoPolarity, PhenoMorphDeviation, RatioNucleousCellSize,Nb,NeighborhoodNames)

    % Heatmap showing OPheno interactions
    figure;
    h = heatmap(CellPhenotypeNames([4,6,7]),CellPhenotypeNames([4,6,7]),CellPhenotypeInteraction([4,6,7],[4,6,7],2),'Colormap',jet,'ColorLimits',[0 1]);
    h.Title = 'Cellular neighborhoods interaction';
    h.XLabel = 'Neighborhoods';
    h.YLabel = 'Neighborhoods'; 
    
    figure; 
    b = bar(categorical(CellPhenotypeNames),PhenoSize); hold on;
    hold on;
    ylabel('Cell phenotype radius (in pixels)');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(F-1,:);
    b.CData(2,:) = colors(F,:);   
    for i=1:F-2        
        b.CData(i+2,:) = colors(i,:);
    end

    figure; 
    b = bar(categorical(CellPhenotypeNames),1-PhenoEccentricity); hold on;
    hold on;
    ylabel('Cell phenotype eccentricity');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(F-1,:);
    b.CData(2,:) = colors(F,:);   
    for i=1:F-2        
        b.CData(i+2,:) = colors(i,:);
    end
    
    figure; 
    b = bar(categorical(CellPhenotypeNames),PhenoPolarity); hold on;
    hold on;
    ylabel('Cell phenotype polarity');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(F-1,:);
    b.CData(2,:) = colors(F,:);   
    for i=1:F-2        
        b.CData(i+2,:) = colors(i,:);
    end
    
    figure; 
    b = bar(categorical(CellPhenotypeNames),RatioNucleousCellSize); hold on;
    hold on;
    ylabel('Cell phenotype eccentricity');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(F-1,:);
    b.CData(2,:) = colors(F,:);   
    for i=1:F-2        
        b.CData(i+2,:) = colors(i,:);
    end

    
    figure; 
    b = bar(categorical(CellPhenotypeNames),PhenoMorphDeviation); hold on;
    hold on;
    ylabel('Cell phenotype eccentricity');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    b.CData(1,:) = colors(F-1,:);
    b.CData(2,:) = colors(F,:);   
    for i=1:F-2        
        b.CData(i+2,:) = colors(i,:);
    end


    figure;
    Phenotype_Abundance(Phenotype_Abundance==0) = 0.0;
    Phenotype_Abundance = Phenotype_Abundance*100;
    b = bar(categorical(NeighborhoodNames),Phenotype_Abundance,'stacked'); hold on;
    hold on;
    ylabel('Relative phenotype abundance');
    colors = double(squeeze(label2rgb(1:F)))/255;    
    for i=1:F       
        b(i).FaceColor = 'flat';
        b(i).CData = colors(i,:);
    end
    legend(CellPhenotypeNames)

end