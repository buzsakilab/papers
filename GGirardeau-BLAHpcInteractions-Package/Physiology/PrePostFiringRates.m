function PrePostFR = PrePostFiringRates(session,presleep,postsleep)
%PrePostFiringRates - Calculates firing rates of individual cells in a session in each of the three states, before and after learning.
%
%  USAGE
%  
%    function PrePostFR = PrePostFiringRates(session,presleep,postsleep)
%
%    session		path to session
%    presleep           presleep epoch index numbered
%    postsleep          postsleep epoch index numbered
%
%  OUTPUT
%
%    PrePostFR : matrix [Rat/session/shank/unit/id/preremrate/postremrate/preswsrate/postswsrate/prewakerate/postwakerate/prerunrate/runrate/postrunrate]
%
%  NOTES
%
%  SEE
%
% 2016 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


PrePostFR=[];

cd(session);
xml=session(end-13:end);
SetCurrentSession(xml);

load('/media/Data-01/All-Rats/sessionindexing.mat');
ratsess = ratsessionindex(strcmp(xmlpath,[session '/']),:)

% Get All Spikes
spikes1=GetSpikeTimes('output','numbered');
spikes2=GetSpikeTimes('output','full');
spikes=[spikes2 spikes1(:,2)];
spikes(spikes(:,3)==0,:)=[];
spikes(spikes(:,3)==1,:)=[];

% Remnumber IDs after removing MUA and artifacts
ids=unique(spikes(:,4));
for ii=1:length(ids)
  spikes(spikes(:,4)==ids(ii),4)=ii;
end
idx=unique(spikes(:,2:4),'rows'); % shank / cell/ ID

if min(idx(:,3))==1 & sum(diff(idx(:,3))-1)==0
  disp('Renumbering OK')
else
  error('Problem with index : aborting')
end

% Brain states
load('States.mat');
remintervals=Rem;
swsintervals=sws;
explointervals=wake;
drowsyintervals=drowsy;

load('runintervals.mat');
if size(runintervals,1)==3
  prerun=runintervals(1,:);
  run=runintervals(2,:);
  postrun=runintervals(3,:);
else
  error('less than 3 run intervals : manual check needed')
end

allPRE=RunIntervals(presleep);
allPOST=RunIntervals(postsleep);
if size(allPRE,1)>1
  allPRE=[allPRE(1,1) allPRE(end,2)];
end
if size(allPOST,1)>1
  allPOST=[allPOST(1,1) allPOST(end,2)];
end

prerem=Restrict(Rem,allPRE);
postrem=Restrict(Rem,allPOST);
presws=Restrict(sws,allPRE);
postsws=Restrict(sws,allPOST);
prewake=Restrict(wake,allPRE);
postwake=Restrict(wake,allPOST);

preremtime=sum(prerem(:,2)-prerem(:,1));
postremtime=sum(postrem(:,2)-postrem(:,1));
preswstime=sum(presws(:,2)-presws(:,1));
postswstime=sum(postsws(:,2)-postsws(:,1));
prewaketime=sum(prewake(:,2)-prewake(:,1));
postwaketime=sum(postwake(:,2)-postwake(:,1));
runtime=sum(run(:,2)-run(:,1));
preruntime=sum(prerun(:,2)-prerun(:,1));
postruntime=sum(postrun(:,2)-postrun(:,1));

allspiketimes=spikes(:,1);
[is.prerun,~,~]=InIntervals(allspiketimes,prerun);
[is.postrun,~,~]=InIntervals(allspiketimes,postrun);
[is.run,~,~]=InIntervals(allspiketimes,run);
[is.prerem,~,~]=InIntervals(allspiketimes,prerem);
[is.presws,~,~]=InIntervals(allspiketimes,presws);
[is.prewake,~,~]=InIntervals(allspiketimes,prewake);
[is.postrem,~,~]=InIntervals(allspiketimes,postrem);
[is.postsws,~,~]=InIntervals(allspiketimes,postsws);
[is.postwake,~,~]=InIntervals(allspiketimes,postwake);

for id=1:size(idx,1)    
  preremrate=sum(spikes(:,4)==id&is.prerem)./preremtime;
  postremrate=sum(spikes(:,4)==id&is.postrem)./postremtime;
  preswsrate=sum(spikes(:,4)==id&is.presws)./preswstime;
  postswsrate=sum(spikes(:,4)==id&is.postsws)./postswstime;
  prewakerate=sum(spikes(:,4)==id&is.prewake)./prewaketime;
  postwakerate=sum(spikes(:,4)==id&is.postwake)./postwaketime;
  prerunrate=sum(spikes(:,4)==id&is.prerun)./preruntime;
  postrunrate=sum(spikes(:,4)==id&is.postrun)./postruntime;
  runrate=sum(spikes(:,4)==id&is.run)./postruntime;
%    fill in the output
  PrePostFR(id,:)=[ratsess idx(id,:) preremrate postremrate preswsrate postswsrate prewakerate postwakerate prerunrate runrate postrunrate];
end

save([xml '-PrePostFR.mat'],'PrePostFR');