function [Neighb_RelAb] = Differences_neighborhood_Ab(Neighborhoods,iImages,Neighb_RelAb)

Neighb_RelAb(iImages,:) = histcounts(Neighborhoods)./sum(sum(histcounts(Neighborhoods)));

figure; boxplot([Neighb_RelAb(:,1),Neighb_RelAb(:,2),Neighb_RelAb(:,3),...
                 Neighb_RelAb(:,4),Neighb_RelAb(:,5),Neighb_RelAb(:,6)],...
                 'Labels',{'Nb1','Nb2','Nb3','Nb4','Nb5','Nb6'}); 


end