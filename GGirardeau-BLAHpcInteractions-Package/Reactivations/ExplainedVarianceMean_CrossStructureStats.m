function ExplainedVarianceMean_CrossStructureStats

% ExplainedVarianceMean_CrossStructureStats - Statitics for EV/REV comparisons across structures.

[censess,cenEVREV,cenp,censtats]=ExplainedVarianceMean_Ultimate('CeCM','EVtype','cross');
[blasess,blaEVREV,blap,blastats]=ExplainedVarianceMean_Ultimate('BLA','EVtype','cross');
[pirsess,pirEVREV,pirp,pirstats]=ExplainedVarianceMean_Ultimate('Pir','EVtype','cross');

%% sws all
bladiff.sws.all=blaEVREV.sws.all.Rall(:,1)-blaEVREV.sws.all.Rall(:,2);
pirdiff.sws.all=pirEVREV.sws.all.Rall(:,1)-pirEVREV.sws.all.Rall(:,2);
cendiff.sws.all=cenEVREV.sws.all.Rall(:,1)-cenEVREV.sws.all.Rall(:,2);

[p,h,st]=kruskalwallis([bladiff.sws.all;pirdiff.sws.all;cendiff.sws.all],[ones(length(bladiff.sws.all),1);ones(length(pirdiff.sws.all),1)*2;ones(length(cendiff.sws.all),1)*3])

bladiff.sws.pyr=blaEVREV.sws.pyr.Rall(:,1)-blaEVREV.sws.pyr.Rall(:,2);
pirdiff.sws.pyr=pirEVREV.sws.pyr.Rall(:,1)-pirEVREV.sws.pyr.Rall(:,2);
cendiff.sws.all=cenEVREV.sws.all.Rall(:,1)-cenEVREV.sws.all.Rall(:,2);
[p,h,st]=kruskalwallis([bladiff.sws.pyr;pirdiff.sws.pyr;cendiff.sws.all],[ones(length(bladiff.sws.pyr),1);ones(length(pirdiff.sws.pyr),1)*2;ones(length(cendiff.sws.all),1)*3])

[p,h,st]=ranksum(bladiff.sws.pyr,pirdiff.sws.pyr)

figure;
plot([ones(length(bladiff.sws.pyr),1);ones(length(pirdiff.sws.pyr),1)*2;ones(length(cendiff.sws.all),1)*3],[bladiff.sws.pyr;pirdiff.sws.pyr;cendiff.sws.all],'k.')