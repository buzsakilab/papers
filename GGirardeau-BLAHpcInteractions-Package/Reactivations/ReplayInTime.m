function [HPtrains,STtrains,R,meanR,tb] = ReplayInTime (session,structure,pre,post,varargin)

%ReplayInTime - Calculates reactivation strength in time during pre and post-sleep epochs for each session.
%
%  USAGE
%
%    [[HPtrains,STtrains,R,meanR,tb] = ReplayInTime (session,structure, pre, post, <options>)
%
%    session		path to session
%    structure          structure name (ex : 'BLA')
%    pre                subsession number for presleep
%    post               subsession number for postsleep
%    <options>          list of property-value pairs, see below
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'binsize'	        binsize for spike trains
%     'savefig'     	save figure : 'on','off'(default)
%     'savevar'		save variable
%     'rzsc'		zscore reactivation strenght (default : 'on')
%     'window'		window for peri-ripple reactivation strength [-2 2]
%     'ctype'	        celltype ('pyr','all')
%    =========================================================================
%
%  OUTPUT
%
%    HPtrains           hippocampal binned spike trains
%    BLAtrains          BLA binned spike trains
%    R                  Reactivation strength
%    meanR              Mean reactivation strength around ripples  
%    tb                 timebins for peri-ripple Rs
%
%  NOTE
%
%  SEE
%
%    See also : ReplayInTime_All, ReplayInTime_Plot, ReplayInTime_LapTypes
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
binsize = 0.05;
savefig = 'off';
savevar = 'off';
rzsc = 'on';
window=[-2 2];
ctype='pyr';

% Check number of inputs
if nargin < 4,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Check varargin
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters  (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'type',
      type = lower(varargin{i+1});
    case 'binsize',
      binsize = (varargin{i+1});
    case 'savevar',
      savevar = (varargin{i+1});
    case 'savefig',
      savefig = (varargin{i+1});
    case 'ripplemod',
      ripplemod = (varargin{i+1});
    case 'zsc'
      rzsc = (varargin{i+1});
    case 'twindow'
      window = (varargin{i+1});
    case 'ctype'
      ctype = (varargin{i+1});      
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
   end
end

load('/media/Data-01/All-Rats/sessionindexing.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
ratsess=ratsessionindex(strcmp([session '/'],xmlpath),1:2); %[rat session]

cd(session);
xml=session(end-13:end);
SetCurrentSession('filename',[session '/' xml]);

%  load([xml '-PoissonRippleMod.mat']);
load('States.mat');
load('runintervals.mat');
load(['/media/Data-01/All-Rats/Structures/' structure '.mat']);
load(['/media/Data-01/All-Rats/Structures/Hpc.mat']);
struc=eval(structure);

% Get pre/post sleep intervals
presleep=RunIntervals(pre);
postsleep=RunIntervals(post);
run=runintervals(2,:);

endlimit=runintervals(end,end);
limits=[0 endlimit];

rip=GetRippleEvents;
rip=rip(:,2);
sleeprip=Restrict(rip,sws);

presws=Restrict(sws,presleep);
postsws=Restrict(sws,postsleep);

% Get hippocampal cells (pyr/all) + index
HPshanks=Hpc(ismember(Hpc(:,1:2),ratsess,'rows'),3);
HPcells=finalType(ismember(finalType(:,1:2),ratsess,'rows')&ismember(finalType(:,3),HPshanks),3:4);
HPspikes=GetSpikeTimes(HPcells,'output','total');
HPindex=unique(HPspikes(:,2:4),'rows');
HPtype=finalType(ismember(finalType(:,1:2),ratsess,'rows')&ismember(finalType(:,3),HPshanks),5); % cell type 1/2
HPispyr=HPtype==1;

% Get all structure spikes + index + finaltype and modulation
STshanks=struc(ismember(struc(:,1:2),ratsess,'rows'),3);
STcells=finalType(ismember(finalType(:,1:2),ratsess,'rows')&ismember(finalType(:,3),STshanks),3:4); %finaltype jsut used as index here
STtype=finalType(ismember(finalType(:,1:2),ratsess,'rows')&ismember(finalType(:,3),STshanks),5);% FinalType for structure cells
STispyr=STtype==1;

if ~isempty(STcells)
  STspikes=GetSpikeTimes(STcells,'output','total');
  STindex=unique(STspikes(:,2:4),'rows');
else
  warning('No cells in target structure : exit function')
  return
end

% Remove extraspikes for session with rec after run-post (ex : Rat08-20130710). Spikes are matrices spiketime/shank/cellmumber/ID
HPspikes(HPspikes(:,1)>=limits(2),:)=[];
STspikes(STspikes(:,1)>=limits(2),:)=[];

% Bin spikes
[HPtrains.raw,bins]=SpikeTrain(HPspikes(:,[1 4]),binsize,limits);
[STtrains.raw,bins]=SpikeTrain(STspikes(:,[1 4]),binsize,limits);

%  Bins.presleep=Restrict(bins,presleep);
%  Bins.run=Restrict(bins,run);
%  Bins.postsleep=Restrict(bins,postsleep);
Bins.presleep=Restrict(bins,presws);
Bins.run=Restrict(bins,run);
Bins.postsleep=Restrict(bins,postsws);

%Construct z-scored Spiketrain matrices for different epochs, pyr only
  
    is.run=InIntervals(bins,runintervals(2,:));
    %  is.presleep=InIntervals(bins,presleep);
    %  is.postsleep=InIntervals(bins,postsleep);
    is.presws=InIntervals(bins,presws);
    is.postsws=InIntervals(bins,postsws);

if strcmp(ctype,'pyr')
    HPtrains.run=zscore(HPtrains.raw(is.run,HPispyr),0,1);
    %  HPtrains.presleep=zscore(HPtrains.raw(is.presleep,HPispyr),0,1);
    %  HPtrains.postsleep=zscore(HPtrains.raw(is.postsleep,HPispyr),0,1);
    HPtrains.presws=zscore(HPtrains.raw(is.presws,HPispyr),0,1);
    HPtrains.postsws=zscore(HPtrains.raw(is.postsws,HPispyr),0,1);

    STtrains.run=zscore(STtrains.raw(is.run,STispyr),0,1);
    %  STtrains.presleep=zscore(STtrains.raw(is.presleep,STispyr),0,1);
    %  STtrains.postsleep=zscore(STtrains.raw(is.postsleep,STispyr),0,1);
    STtrains.presws=zscore(STtrains.raw(is.presws,STispyr),0,1);
    STtrains.postsws=zscore(STtrains.raw(is.postsws,STispyr),0,1);
elseif strcmp(ctype,'all')
    HPtrains.run=zscore(HPtrains.raw(is.run,:),0,1);
    %  HPtrains.presleep=zscore(HPtrains.raw(is.presleep,HPispyr),0,1);
    %  HPtrains.postsleep=zscore(HPtrains.raw(is.postsleep,HPispyr),0,1);
    HPtrains.presws=zscore(HPtrains.raw(is.presws,:),0,1);
    HPtrains.postsws=zscore(HPtrains.raw(is.postsws,:),0,1);

    STtrains.run=zscore(STtrains.raw(is.run,:),0,1);
    %  STtrains.presleep=zscore(STtrains.raw(is.presleep,STispyr),0,1);
    %  STtrains.postsleep=zscore(STtrains.raw(is.postsleep,STispyr),0,1);
    STtrains.presws=zscore(STtrains.raw(is.presws,:),0,1);
    STtrains.postsws=zscore(STtrains.raw(is.postsws,:),0,1);
end

%RUN z-scored correlation matrix
corrM.bla=STtrains.run'*STtrains.run/sum(is.run);
corrM.hpc=((1/sum(is.run))*HPtrains.run'*HPtrains.run);
corrM.hpbla=(1/sum(is.run))*HPtrains.run'*STtrains.run;
corrM.giant=(1/sum(is.run))*[HPtrains.run STtrains.run]'*[HPtrains.run STtrains.run];

f1=figure;
subplot(1,3,1)
imagesc(corrM.bla);
subplot(1,3,2)
imagesc(corrM.hpc);
subplot(1,3,3)
imagesc(corrM.hpbla);

%  figure;
%  imagesc(corrM.giant);

%WAke (raw) matrix reactivations during sleep
%  for i=1:sum(is.presleep)
%    R.bla.presleep(i)=STtrains.presleep(i,:)*corrM.bla*STtrains.presleep(i,:)'-STtrains.presleep(i,:)*STtrains.presleep(i,:)';
%    R.hpc.presleep(i)=HPtrains.presleep(i,:)*corrM.hpc*HPtrains.presleep(i,:)'-HPtrains.presleep(i,:)*HPtrains.presleep(i,:)';
%    R.cross.presleep(i)=HPtrains.presleep(i,:)*corrM.hpbla*STtrains.presleep(i,:)';
%  end
%  for i=1:sum(is.postsleep)
%    R.bla.postsleep(i)=STtrains.postsleep(i,:)*corrM.bla*STtrains.postsleep(i,:)'-STtrains.postsleep(i,:)*STtrains.postsleep(i,:)';
%    R.hpc.postsleep(i)=HPtrains.postsleep(i,:)*corrM.hpc*HPtrains.postsleep(i,:)'-HPtrains.postsleep(i,:)*HPtrains.postsleep(i,:)';
%    R.cross.postsleep(i)=HPtrains.postsleep(i,:)*corrM.hpbla*STtrains.postsleep(i,:)';
%  end
for i=1:sum(is.presws)
  R.bla.presws(i)=STtrains.presws(i,:)*corrM.bla*STtrains.presws(i,:)'-STtrains.presws(i,:)*STtrains.presws(i,:)';
  R.hpc.presws(i)=HPtrains.presws(i,:)*corrM.hpc*HPtrains.presws(i,:)'-HPtrains.presws(i,:)*HPtrains.presws(i,:)';
  R.cross.presws(i)=HPtrains.presws(i,:)*corrM.hpbla*STtrains.presws(i,:)';
end
for i=1:sum(is.postsws)
  R.bla.postsws(i)=STtrains.postsws(i,:)*corrM.bla*STtrains.postsws(i,:)'-STtrains.postsws(i,:)*STtrains.postsws(i,:)';
  R.hpc.postsws(i)=HPtrains.postsws(i,:)*corrM.hpc*HPtrains.postsws(i,:)'-HPtrains.postsws(i,:)*HPtrains.postsws(i,:)';
  R.cross.postsws(i)=HPtrains.postsws(i,:)*corrM.hpbla*STtrains.postsws(i,:)';
end

%%%%%%%%%%%% Optional Zscore
if strcmp(rzsc,'on')
%    R.bla.presleep=zscore(R.bla.presleep);
%    R.hpc.presleep=zscore(R.hpc.presleep);
%    R.cross.presleep=zscore(R.cross.presleep);
%    R.bla.postsleep=zscore(R.bla.postsleep);
%    R.hpc.postsleep=zscore(R.hpc.postsleep);
%    R.cross.postsleep=zscore(R.cross.postsleep);
  R.bla.presleep=zscore(R.bla.presws);
  R.hpc.presleep=zscore(R.hpc.presws);
  R.cross.presleep=zscore(R.cross.presws);
  R.bla.postsleep=zscore(R.bla.postsws);
  R.hpc.postsleep=zscore(R.hpc.postsws);
  R.cross.postsleep=zscore(R.cross.postsws);
end

presleeprip=Restrict(sleeprip,presleep);
postsleeprip=Restrict(sleeprip,postsleep);

%  % Plotting
%  f2=figure;
%  subplot(2,1,1)
%  plot(Bins.presleep,R.bla.presleep,'k');hold on;
%  PlotHVLines(presleeprip,'v','r');
%  xlabel('Cross - Presleep');
%  subplot(2,1,2)
%  plot(Bins.postsleep,R.bla.postsleep,'k');hold on;
%  PlotHVLines(postsleeprip,'v','r');
%  xlabel('Cross - Presleep');
%  suptitle(xml);

% Pre/Post Peri-ripple no PCA reactivations
[syncR.bla.pre,idc.bla.pre]=Sync([Bins.presleep R.bla.presleep'],presleeprip,'durations',window);
[syncR.bla.post,idc.bla.post]=Sync([Bins.postsleep R.bla.postsleep'],postsleeprip,'durations',window);
[syncR.hpc.pre,idc.hpc.pre]=Sync([Bins.presleep R.hpc.presleep'],presleeprip,'durations',window);
[syncR.hpc.post,idc.hpc.post]=Sync([Bins.postsleep R.hpc.postsleep'],postsleeprip,'durations',window);
[syncR.cross.pre,idc.cross.pre]=Sync([Bins.presleep R.cross.presleep'],presleeprip,'durations',window);
[syncR.cross.post,idc.cross.post]=Sync([Bins.postsleep R.cross.postsleep'],postsleeprip,'durations',window);

[meanR.bla.pre,err.bla.pre,tb]=SyncHist(syncR.bla.pre,idc.bla.pre,'mode','mean','durations',window,'smooth',2);
[meanR.bla.post,err.bla.post,tb]=SyncHist(syncR.bla.post,idc.bla.post,'mode','mean','durations',window,'smooth',2);
[meanR.hpc.pre,err.hpc.pre,tb]=SyncHist(syncR.hpc.pre,idc.hpc.pre,'mode','mean','durations',window,'smooth',2);
[meanR.hpc.post,err.hpc.post,tb]=SyncHist(syncR.hpc.post,idc.hpc.post,'mode','mean','durations',window,'smooth',2);
[meanR.cross.pre,err.cross.pre,tb]=SyncHist(syncR.cross.pre,idc.cross.pre,'mode','mean','durations',window,'smooth',2,'error','sem');
[meanR.cross.post,err.cross.post,tb]=SyncHist(syncR.cross.post,idc.cross.post,'mode','mean','durations',window,'smooth',2,'error','sem');

f2=figure('Position',[605 492 1151 342]);hold on;
subplot(1,3,1);hold on;
plot(tb,meanR.bla.pre,'b');
plot(tb,meanR.bla.post,'r');
xlabel('BLA');
subplot(1,3,2);hold on;
plot(tb,meanR.hpc.pre,'b');
plot(tb,meanR.hpc.post,'r');
xlabel('Hpc');
subplot(1,3,3);hold on;
plot(tb,meanR.cross.pre,'b');
plot(tb,meanR.cross.post,'r');
xlabel('Cross');
suptitle(xml);

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  %%%% PCA
%  [coeffs.bla,latent.bla,expl.bla]=pcacov(corrM.bla);
%  [coeffs.hpc,latent.hpc,expl.hpc]=pcacov(corrM.hpc);
%  [coeffs.giant,latent.giant,expl.giant]=pcacov(corrM.giant);
%  
%  %  figure;
%  %  subplot(1,3,1);
%  %  stem(coeffs.bla(:,1));
%  %  xlabel('BLA')
%  %  subplot(1,3,2);
%  %  stem(coeffs.hpc(:,1));
%  %  xlabel('hpc');
%  %  subplot(1,3,3);
%  %  stem(coeffs.giant(:,1));
%  %  xlabel('cross');
%  
%  % Intrastructure Reactivations per PC
%  for i=1:size(coeffs.bla,2)
%    Sc.bla.presleep(:,i)=STtrains.presleep*coeffs.bla(:,i);
%    Rpc.bla.presleep(:,i)=Sc.bla.presleep(:,i).^2-STtrains.presleep.^2*coeffs.bla(:,i).^2;
%    Sc.bla.postsleep(:,i)=STtrains.postsleep*coeffs.bla(:,i);
%    Rpc.bla.postsleep(:,i)=Sc.bla.postsleep(:,i).^2-STtrains.postsleep.^2*coeffs.bla(:,i).^2;
%  end
%  
%  for i=1:size(coeffs.hpc,2)
%    Sc.hpc.presleep(:,i)=HPtrains.presleep*coeffs.hpc(:,i);
%    Rpc.hpc.presleep(:,i)=Sc.hpc.presleep(:,i).^2-HPtrains.presleep.^2*coeffs.hpc(:,i).^2;
%    Sc.hpc.postsleep(:,i)=HPtrains.postsleep*coeffs.hpc(:,i);
%    Rpc.hpc.postsleep(:,i)=Sc.hpc.postsleep(:,i).^2-HPtrains.postsleep.^2*coeffs.hpc(:,i).^2; 
%  end
%  
%  for i=1:size(coeffs.giant,2)
%    Sc.giant.presleep(:,i)=[HPtrains.presleep STtrains.presleep]*coeffs.giant(:,i);
%    Rpc.giant.presleep(:,i)=Sc.giant.presleep(:,i).^2-[HPtrains.presleep STtrains.presleep].^2*coeffs.giant(:,i).^2;
%    Sc.giant.postsleep(:,i)=[HPtrains.postsleep STtrains.postsleep]*coeffs.giant(:,i);
%    Rpc.giant.postsleep(:,i)=Sc.giant.postsleep(:,i).^2-[HPtrains.postsleep STtrains.postsleep].^2*coeffs.giant(:,i).^2; 
%  end
%  
%  %z-score optionnel des reactivations
%  if strcmp(rzsc,'on')
%    Rpc.bla.presleep=zscore(Rpc.bla.presleep,0,1);
%    Rpc.bla.postsleep=zscore(Rpc.bla.postsleep,0,1);
%    Rpc.hpc.presleep=zscore(Rpc.hpc.presleep,0,1);
%    Rpc.hpc.postsleep=zscore(Rpc.hpc.postsleep,0,1);
%    Rpc.giant.presleep=zscore(Rpc.giant.presleep,0,1);
%    Rpc.giant.postsleep=zscore(Rpc.giant.postsleep,0,1);
%  end
%    
%  % Signal vs noise PC
%  vMax.bla=(1+sqrt(size(STtrains.postsleep,2)/size(STtrains.postsleep,1)))^2;
%  pcIndex.bla=find(latent.bla>vMax.bla);
%  vMax.hpc=(1+sqrt(size(HPtrains.postsleep,2)/size(HPtrains.postsleep,1)))^2;
%  pcIndex.hpc=find(latent.hpc>vMax.hpc);
%  vMax.giant=(1+sqrt(size([HPtrains.postsleep STtrains.postsleep],2)/size([HPtrains.postsleep STtrains.postsleep],1)))^2;
%  pcIndex.giant=find(latent.giant>vMax.giant);
%  
%  % Peri-ripple reactivations for PC 1 to 3  for intra-BLA and intra-HPcells
%  for i=1:3
%    [syncRpc.bla.pre{i},idc.bla.pre(:,i)]=Sync([Bins.presleep Rpc.bla.presleep(:,i)],presleeprip,'durations',window);
%    [syncRpc.bla.post{i},idc.bla.post(:,i)]=Sync([Bins.postsleep Rpc.bla.postsleep(:,i)],postsleeprip,'durations',window);
%    [meanRpc.bla.pre(:,i),err.bla.pre(:,i),tb]=SyncHist(syncRpc.bla.pre{i},idc.bla.pre(:,i),'mode','mean','durations',window,'smooth',2);
%    [meanRpc.bla.post(:,i),err.bla.post(:,i),tb]=SyncHist(syncRpc.bla.post{i},idc.bla.post(:,i),'mode','mean','durations',window,'smooth',2);
%  
%    [syncRpc.hpc.pre{i},idc.hpc.pre(:,i)]=Sync([Bins.presleep Rpc.hpc.presleep(:,i)],presleeprip,'durations',window);
%    [syncRpc.hpc.post{i},idc.hpc.post(:,i)]=Sync([Bins.postsleep Rpc.hpc.postsleep(:,i)],postsleeprip,'durations',window);
%    [meanRpc.hpc.pre(:,i),err.hpc.pre(:,i),tb]=SyncHist(syncRpc.hpc.pre{i},idc.hpc.pre(:,i),'mode','mean','durations',window,'smooth',2);
%    [meanRpc.hpc.post(:,i),err.hpc.post(:,i),tb]=SyncHist(syncRpc.hpc.post{i},idc.hpc.post(:,i),'mode','mean','durations',window,'smooth',2);
%    
%    [syncRpc.giant.pre{i},idc.giant.pre(:,i)]=Sync([Bins.presleep Rpc.giant.presleep(:,i)],presleeprip,'durations',window);
%    [syncRpc.giant.post{i},idc.giant.post(:,i)]=Sync([Bins.postsleep Rpc.giant.postsleep(:,i)],postsleeprip,'durations',window);
%    [meanRpc.giant.pre(:,i),err.giant.pre(:,i),tb]=SyncHist(syncRpc.giant.pre{i},idc.giant.pre(:,i),'mode','mean','durations',window,'smooth',2);
%    [meanRpc.giant.post(:,i),err.giant.post(:,i),tb]=SyncHist(syncRpc.giant.post{i},idc.giant.post(:,i),'mode','mean','durations',window,'smooth',2);
%  end
%  
%  f3=figure('Position',[703   137   927   816]);hold on;
%  subplot(3,3,1);hold on;
%  plot(tb,meanRpc.bla.pre(:,1),'b');
%  plot(tb,meanRpc.bla.post(:,1),'r');
%  ylabel('Reactivations PC1');
%  subplot(3,3,2);hold on;
%  plot(tb,meanRpc.hpc.pre(:,1),'b');
%  plot(tb,meanRpc.hpc.post(:,1),'r');
%  subplot(3,3,3);hold on;
%  plot(tb,meanRpc.giant.pre(:,1),'b');
%  plot(tb,meanRpc.giant.post(:,1),'r');
%  subplot(3,3,4);hold on;
%  plot(tb,meanRpc.bla.pre(:,2),'b');
%  plot(tb,meanRpc.bla.post(:,2),'r');
%  ylabel('Reactivations PC2');
%  subplot(3,3,5);hold on;
%  plot(tb,meanRpc.hpc.pre(:,2),'b');
%  plot(tb,meanRpc.hpc.post(:,2),'r');
%  subplot(3,3,6);hold on;
%  plot(tb,meanRpc.giant.pre(:,2),'b');
%  plot(tb,meanRpc.giant.post(:,2),'r'); 
%  subplot(3,3,7);hold on;
%  plot(tb,meanRpc.bla.pre(:,3),'b');
%  plot(tb,meanRpc.bla.post(:,3),'r');
%  xlabel('BLA or Struc');
%  ylabel('Reactivations PC3')
%  subplot(3,3,8);hold on;
%  plot(tb,meanRpc.hpc.pre(:,3),'b');
%  plot(tb,meanRpc.hpc.post(:,3),'r');
%  xlabel('Hpc');
%  subplot(3,3,9);hold on;
%  plot(tb,meanRpc.giant.pre(:,3),'b');
%  plot(tb,meanRpc.giant.post(:,3),'r'); 
%  xlabel('Hpc and BLA/struc giant');
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % Cross structure reactivations per PC, partial fit.
%  for j=1:size(coeffs.giant,2)
%    proj=[];
%    proj=coeffs.giant(:,j)*coeffs.giant(:,j)';
%    proj=proj(sum(HPispyr)+1:end,1:sum(HPispyr));
%    for i=1:sum(is.presleep)
%      Rpc.cross.pre(j,i)=STtrains.presleep(i,:)*proj*HPtrains.presleep(i,:)';
%    end
%    for i=1:sum(is.postsleep)
%      Rpc.cross.post(j,i)=STtrains.postsleep(i,:)*proj*HPtrains.postsleep(i,:)';
%    end
%  end
%  
%  if strcmp(rzsc,'on')
%    Rpc.cross.pre=zscore(Rpc.cross.pre,0,2);
%    Rpc.cross.post=zscore(Rpc.cross.post,0,2);
%  end
%  
%  %  figure;suptitle('CrossPC')
%  %  for i=1:3
%  %    subplot(3,1,i);hold on;
%  %    plot(Bins.presleep,Rpc.cross.pre(i,:),'b');
%  %    plot(Bins.postsleep,Rpc.cross.post(i,:),'r');
%  %  end
%  
%  % Peri-ripple reactivations for PC 1 to 3  for intra-BLA and intra-HPcells
%  for i=1:9
%    [syncRpc.cross.pre{i},idc.cross.pre(:,i)]=Sync([Bins.presleep Rpc.cross.pre(i,:)'],presleeprip,'durations',window);
%    [syncRpc.cross.post{i},idc.cross.post(:,i)]=Sync([Bins.postsleep Rpc.cross.post(i,:)'],postsleeprip,'durations',window);
%    [meanRpc.cross.pre(:,i),err.cross.pre(:,i),tb]=SyncHist(syncRpc.cross.pre{i},idc.cross.pre(:,i),'mode','mean','durations',window,'smooth',2);
%    [meanRpc.cross.post(:,i),err.cross.post(:,i),tb]=SyncHist(syncRpc.cross.post{i},idc.cross.post(:,i),'mode','mean','durations',window,'smooth',2);
%  end
%  
%  %  figure;hold on;
%  %  subplot(3,1,1);hold on;
%  %  plot(tb,meanRpc.cross.pre(:,1),'b');
%  %  plot(tb,meanRpc.cross.post(:,1),'r');
%  %  subplot(3,1,2);hold on;
%  %  plot(tb,meanRpc.cross.pre(:,2),'b');
%  %  plot(tb,meanRpc.cross.post(:,2),'r');
%  %  subplot(3,1,3);hold on;
%  %  plot(tb,meanRpc.cross.pre(:,3),'b');
%  %  plot(tb,meanRpc.cross.post(:,3),'r');
%  
%  f4=figure('Position',[703 137 927 816]);hold on;
%  for i=1:9
%    SquareSubplot(9,i);hold on;
%    plot(tb,meanRpc.cross.pre(:,i),'b');
%    plot(tb,meanRpc.cross.post(:,i),'r');
%    xlabel(['PC ' num2str(i)])
%  end
%  suptitle('PCA-CrossStructure-PartialFit-PeriSPWR');
%  
if strcmp(savevar,'on')
%    save([xml '-ReplayInTime-' structure '-binsize' num2str(binsize) '-zsc' rzsc '.mat'],'HPtrains','STtrains','R','Rpc','meanR','meanRpc','tb','corrM','bins','Bins')
  save([xml '-ReplayInTime-' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.mat'],'HPtrains','STtrains','R','meanR','tb','corrM','bins','Bins')
end

if strcmp(savefig,'on')
  cd('/media/Data-01/All-Rats/AllRats-ReplayInTime');
  plot2svg(['TimeReplay' xml '-CorrMatrices-' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.svg'],f1); %%
  plot2svg(['TimeReplay' xml '-RipplePrePost-NoPCA-' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.svg'],f2);
%    plot2svg(['TimeReplay' xml '-RipplePrePost-ClassicPCA-PC1a3-' structure '-binsize' num2str(binsize) '-zsc' rzsc '.svg'],f3);
%    plot2svg(['TimeReplay' xml '-RipplePrePost-CrossPCA-PC1a9-' structure '-binsize' num2str(binsize) '-zsc' rzsc '.svg'],f4);
  close all
end
