function [] = Parameters_Visualization_Phenotypes(CellPhenotypeNames,PhenoSize,F,...
                                                Phenotype_Abundance,PhenoEccentricity,Nb,NeighborhoodNames)

    figure; 
    b = bar(categorical(CellPhenotypeNames),PhenoSize); hold on;
    hold on;
    ylabel('Cell phenotype radius (in pixels)');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    for i=1:F        
        b.CData = colors(i,:);
    end

    figure; 
    b = bar(categorical(CellPhenotypeNames),1-PhenoEccentricity); hold on;
    hold on;
    ylabel('Cell phenotype eccentricity');
    colors = double(squeeze(label2rgb(1:F)))/255;
    b.FaceColor = 'flat';
    for i=1:F        
        b.CData(i,:) = colors(i,:);
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