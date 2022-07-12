function [step] = Modeling_Neighborhoods(BasefolderName,LoadfolderName,...
                           sizeImages,subject_type,subject_number,verbose)
% Variable initialization
Perc_Pheno=[];
RandIndicesPrev=[];
saveOptimizationLossNeigh = [0];
saveOptimizationNumberOfPixelsNeigh = [0];
load(['Synthetic_microenvironment/parameters/CellularNeighborhoods_',subject_type,'_',num2str(subject_number),'.mat']);
Neighborhoods_Interaction = Neighborhoods_Interaction+0.001

% Normalize Neighborhoods
Neighborhood_Abundances = (Neighborhood_Abundances./sum(Neighborhood_Abundances)); 
Neighborhood_Abundances(Neighborhood_Abundances==0) = .0001;

% Initialize neighborhoods mask ground truth
GT_Nb = reshape(randsample(Nb,sizeImages(1)*sizeImages(2),true,ones(Nb,1)'),sizeImages);

% Initialize list of painted pixels
excludedIndices = [];

% Contextual Limitation to assign Neighborhoods.
Contxt_limit = max(NeighborhoodMinSize);

% Obtain indices that are inside the limits of the image.
Indices = Modeling_indices_inside_limits(GT_Nb,3);

% Randomice indices to update pixels in a random order
RandIndices = Indices(randperm(length(Indices)));

% Update a 0.1% of pixels each step
stepLength = 0:length(RandIndices)*0.001:length(RandIndices); 
step = 0;

% Size of the ground-truth
[Imr,Imc,Imz]=size(GT_Nb);
 
% Start optimization
fprintf(['.']); lastIndex=0;iteration=1;
while ~isempty(RandIndices)
          
    % Display the percentages of Communities(% of total) and Phenotypes(% of community). 
    GT_Nb_ctx=GT_Nb(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit);
    GT_Nb_Ab = histcounts([GT_Nb_ctx(:)',[1:Nb]],'BinMethod', 'Integers')';            
    GT_Nb_Ab = GT_Nb_Ab./sum(GT_Nb_Ab(:))*100;
    Diff_Nb = GT_Nb_Ab' - Neighborhood_Abundances*100;        
    fprintf(['Loss:',num2str(sum(abs(Diff_Nb))),'_Pixels:',num2str(length(RandIndices)),'\n'])                

    % Save optimization information
    saveOptimizationLossNeigh = cat(1,saveOptimizationLossNeigh,sum(abs(Diff_Nb)));
    saveOptimizationNumberOfPixelsNeigh = cat(1,saveOptimizationNumberOfPixelsNeigh,length(RandIndices));
    
    % Eliminate iterations if pixels were already assigned.
    RandIndices=RandIndices(~ismember(RandIndices,excludedIndices));
    RandIndices(randperm(numel(RandIndices))) = RandIndices;

    % In each iteration we reduce the interaction between neighborhoods.
    fprintf([num2str(length(RandIndices)),' '])
    if length(RandIndices)>5000               
        if sum(abs(Diff_Nb))>4                    
            Numstep=[1000];
            Neighborhoods_Interaction = Neighborhoods_Interaction+eye(length(Neighborhood_Abundances))*0.000001;
        else
            Neighborhoods_Interaction = Neighborhoods_Interaction+eye(length(Neighborhood_Abundances))*0.000001;
            Numstep=[1000];
        end
    else
        Numstep=length(RandIndices); % Just use the steps necessary to finish the optimization
    end  
    
    % Check if optimization has stopped
    if (lastIndex==length(RandIndices) & length(RandIndices)<length(Indices)) || iteration>=1500
        optimization_is_blocked=true; 
        iteration=975;
    else
        optimization_is_blocked=false;
    end
    lastIndex=length(RandIndices);
    iteration= iteration+1;
    
    % Iterate to assign Neighborhoods
    for i=1:floor(Numstep)
        if ismembc(RandIndices(i),excludedIndices); continue; end            

        % Obtain row/column and check if it is a boundary
        [r,c] = ind2sub([sizeImages(1),sizeImages(2)],RandIndices(i));                                    

        % Calculate update rule                           
        Neighborhoods_Update = (Neighborhoods_Interaction.*(((Neighborhood_Abundances*100)'./GT_Nb_Ab).^2)');                
        Contxt_limit_Rand = round(Contxt_limit-rand*(Contxt_limit*0.25));
        SelblockBIG = GT_Nb(max(r-Contxt_limit_Rand*1,Contxt_limit_Rand*1): ...
                            min(r+Contxt_limit_Rand*1,sizeImages(1)-Contxt_limit_Rand*1),...
                            max(c-Contxt_limit_Rand*1,Contxt_limit_Rand*1):...
                            min(c+Contxt_limit_Rand*1,sizeImages(2)-Contxt_limit_Rand*1));
        SelblockSMALL = GT_Nb(round(max(r-Contxt_limit_Rand/4,Contxt_limit_Rand/4)):...
                              round(min(r+Contxt_limit_Rand/4,sizeImages(1)-Contxt_limit_Rand/4)),...
                              round(max(c-Contxt_limit_Rand/4,Contxt_limit_Rand/4)):...
                              round(min(c+Contxt_limit_Rand/4,sizeImages(2)-Contxt_limit_Rand/4)));

        for nb_i=1:Nb; SelblockPixls(nb_i)=sum(SelblockSMALL(:)==nb_i); end
        SelblockPixls=SelblockPixls./length(SelblockSMALL(:));  
        SelblockPixls= SelblockPixls+0.001;
        if (max(SelblockPixls)>0.4  & sum(abs(Diff_Nb))<2) || optimization_is_blocked%(mode_rand(SelblockSMALL(:))))<2
            
            % Let's assign the new neighborhood to the pixel with a certain probability.
            winner = mode_rand(SelblockSMALL(:));  

            % Create size of community
            SE = strel('disk',1,0); 
            [Rmask,Cmask]=find(SE.Neighborhood==1);
            Rmask=Rmask+r-1-1; Cmask=Cmask+c-1-1;

            maskInd=Rmask+Imr.*(Cmask-1)+Imr.*Imc.*(Imz-1);                
            GT_Nb(maskInd) = winner;              

            if max(SelblockPixls)>0.95 || length(RandIndices) <0.2*sizeImages(1)*sizeImages(2) || optimization_is_blocked
                % Save indices avoid repetition of painted communities.
                excludedIndices = [excludedIndices; Rmask+Imr.*(Cmask-1)];     
            end
        else                    
          
            if length(RandIndices) <0.2*sizeImages(1)*sizeImages(2)
                winner = mode_rand(SelblockBIG(:)); % Calc Community Mode and Slect Probability.                                 
                GT_Nb(r,c) = winner; % Set the new Anatomical Region of the ARmap
                if length(RandIndices)<0.2*sizeImages(1)*sizeImages(2) || optimization_is_blocked
                    excludedIndices = [excludedIndices; RandIndices(i)];
                end
            else
                probComm = mean(Neighborhoods_Update(SelblockBIG(:),:)); % Calc Community Mode and Slect Probability.         
                [~,winner] = max(probComm);

                GT_Nb(r-1:r+1,c-1:c+1) = winner; % Set the new Anatomical Region of the ARmap
            end                   
        end
    end  
    
    % IF true, save tissue neighborhoods at the actual step.
    if verbose==true
        % Save step by step.
        step = step + 1;       
        if mod(step,50)==1
            Modeling_Neighborhoods_Save(BasefolderName,GT_Nb,subject_type,subject_number,step,saveOptimizationLossNeigh,saveOptimizationNumberOfPixelsNeigh)
        end
    end
    
end
  
% Save final step.
step = step + 1;  
GT_Nb = modefilt(GT_Nb,[NeighborhoodMinSize(1)*3+1, NeighborhoodMinSize(1)*3+1]);
Modeling_Neighborhoods_Save(BasefolderName,GT_Nb,subject_type,subject_number,step,saveOptimizationLossNeigh,saveOptimizationNumberOfPixelsNeigh)                

end