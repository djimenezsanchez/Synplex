function [interactivity_global] = Differences_between_CellNeighb(iImages,Neighborhoods,interactivity_global)

% Get interactions    
Neighborhoods_BW = imdilate(edge(Neighborhoods,'Canny'),strel('disk',4));
Neighborhoods = Neighborhoods.*Neighborhoods_BW;
interactivity_local = zeros(max(Neighborhoods(:)),max(Neighborhoods(:)));
for x=1:20:size(Neighborhoods,1)-20
    for y=1:20:size(Neighborhoods,2)-20
        count = zeros(max(Neighborhoods(:)),1);
        m = Neighborhoods(x:x+20,y:y+20);
        for chan=1:max(Neighborhoods(:))                
            count(chan)= sum(sum(m==chan));
        end
        if sum(count>1)==2
            n= find(count);
            interactivity_local(n(1),n(2)) = interactivity_local(n(1),n(2))+1;
            interactivity_local(n(2),n(1)) = interactivity_local(n(2),n(1))+1;
        end
    end
end
interactivity_local = interactivity_local/sum(interactivity_local(:));
for chan1=1:max(Neighborhoods(:))
    for chan2=1:max(Neighborhoods(:))
        interactivity_global(chan1,chan2,iImages) = interactivity_global(chan1,chan2,iImages)+interactivity_local(chan1,chan2);
    end
end                                            

end