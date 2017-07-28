function [HPtrains,STtrains,R,meanR,tb] = ReplayInTime_LapTypes (session,structure,pre,post,varargin)

%ReplayInTime_LapTypes - Calculates reactivation strength of the airpuff vs the safe trajectory correlation matrices in time during pre and post-sleep epochs for each session.
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

load([xml '-LapType.mat']);

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

Bins.presleep=Restrict(bins,presws);
Bins.run=Restrict(bins,run);
Bins.postsleep=Restrict(bins,postsws);

%Construct z-scored Spiketrain matrices for different epochs, pyr only
is.run.safe=InIntervals(bins,safelaps.run);
is.run.ap=InIntervals(bins,aplaps.run);
is.run.all=InIntervals(bins,runintervals(2,:));
is.presws=InIntervals(bins,presws);
is.postsws=InIntervals(bins,postsws);

if strcmp(ctype,'pyr')
    HPtrains.run.safe=zscore(HPtrains.raw(is.run.safe,HPispyr),0,1);
    HPtrains.run.ap=zscore(HPtrains.raw(is.run.ap,HPispyr),0,1);
    HPtrains.presws=zscore(HPtrains.raw(is.presws,HPispyr),0,1);
    HPtrains.postsws=zscore(HPtrains.raw(is.postsws,HPispyr),0,1);

    STtrains.run.safe=zscore(STtrains.raw(is.run.safe,STispyr),0,1);
    STtrains.run.ap=zscore(STtrains.raw(is.run.ap,STispyr),0,1);
    STtrains.presws=zscore(STtrains.raw(is.presws,STispyr),0,1);
    STtrains.postsws=zscore(STtrains.raw(is.postsws,STispyr),0,1);
elseif strcmp(ctype,'all')
    HPtrains.run.safe=zscore(HPtrains.raw(is.run.safe,:),0,1);
    HPtrains.run.ap=zscore(HPtrains.raw(is.run.ap,:),0,1);
    HPtrains.presws=zscore(HPtrains.raw(is.presws,:),0,1);
    HPtrains.postsws=zscore(HPtrains.raw(is.postsws,:),0,1);

    STtrains.run.safe=zscore(STtrains.raw(is.run.safe,:),0,1);
    STtrains.run.ap=zscore(STtrains.raw(is.run.ap,:),0,1);
    STtrains.presws=zscore(STtrains.raw(is.presws,:),0,1);
    STtrains.postsws=zscore(STtrains.raw(is.postsws,:),0,1);
end

%RUN z-scored correlation matrix
corrM.hpbla.safe=(1/sum(is.run.safe))*HPtrains.run.safe'*STtrains.run.safe;
corrM.hpbla.ap=(1/sum(is.run.ap))*HPtrains.run.ap'*STtrains.run.ap;

f1=figure;
subplot(1,2,1)
imagesc(corrM.hpbla.safe);
xlabel('safe')
subplot(1,2,2)
imagesc(corrM.hpbla.ap);
xlabel('airpuff')

for i=1:sum(is.presws)
  R.cross.presws.safe(i)=HPtrains.presws(i,:)*corrM.hpbla.safe*STtrains.presws(i,:)';
  R.cross.presws.ap(i)=HPtrains.presws(i,:)*corrM.hpbla.ap*STtrains.presws(i,:)';
end
for i=1:sum(is.postsws)
  R.cross.postsws.safe(i)=HPtrains.postsws(i,:)*corrM.hpbla.safe*STtrains.postsws(i,:)';
  R.cross.postsws.ap(i)=HPtrains.postsws(i,:)*corrM.hpbla.ap*STtrains.postsws(i,:)';
end

%%%%%%%%%%%% Optional Zscore
if strcmp(rzsc,'on')
  R.cross.presleep.safe=zscore(R.cross.presws.safe);
  R.cross.presleep.ap=zscore(R.cross.presws.ap);
  R.cross.postsleep.safe=zscore(R.cross.postsws.safe);
  R.cross.postsleep.ap=zscore(R.cross.postsws.ap);
end

presleeprip=Restrict(sleeprip,presleep);
postsleeprip=Restrict(sleeprip,postsleep);

% Pre/Post Peri-ripple no PCA reactivations
[syncR.cross.pre.safe,idc.cross.pre.safe]=Sync([Bins.presleep R.cross.presleep.safe'],presleeprip,'durations',window);
[syncR.cross.pre.ap,idc.cross.pre.ap]=Sync([Bins.presleep R.cross.presleep.ap'],presleeprip,'durations',window);

[syncR.cross.post.safe,idc.cross.post.safe]=Sync([Bins.postsleep R.cross.postsleep.safe'],postsleeprip,'durations',window);
[syncR.cross.post.ap,idc.cross.post.ap]=Sync([Bins.postsleep R.cross.postsleep.ap'],postsleeprip,'durations',window);

[meanR.cross.pre.safe,err.cross.pre.safe,tb]=SyncHist(syncR.cross.pre.safe,idc.cross.pre.safe,'mode','mean','durations',window,'smooth',2,'error','sem');
[meanR.cross.pre.ap,err.cross.pre.ap,tb]=SyncHist(syncR.cross.pre.ap,idc.cross.pre.ap,'mode','mean','durations',window,'smooth',2,'error','sem');

[meanR.cross.post.safe,err.cross.post.safe,tb]=SyncHist(syncR.cross.post.safe,idc.cross.post.safe,'mode','mean','durations',window,'smooth',2,'error','sem');
[meanR.cross.post.ap,err.cross.post.ap,tb]=SyncHist(syncR.cross.post.ap,idc.cross.post.ap,'mode','mean','durations',window,'smooth',2,'error','sem');

f2=figure('Position',[605 492 1151 342]);hold on;
subplot(1,2,1);hold on;
plot(tb,meanR.cross.pre.safe,'b');
plot(tb,meanR.cross.post.safe,'r');
xlabel('SAFE');
subplot(1,2,2);hold on;
plot(tb,meanR.cross.pre.ap,'b');
plot(tb,meanR.cross.post.ap,'r');
xlabel('AP')
suptitle(xml);

%  
if strcmp(savevar,'on')
  save([xml '-ReplayInTime-LapTypes-' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.mat'],'HPtrains','STtrains','R','meanR','tb','corrM','bins','Bins')
end

if strcmp(savefig,'on')
  cd('/media/Data-01/All-Rats/AllRats-ReplayInTime');
  plot2svg(['TimeReplay' xml '-CorrMatrices-LapTypes' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.svg'],f1); %%
  plot2svg(['TimeReplay' xml '-RipplePrePost-NoPCA-LapTypes' structure '-binsize' num2str(binsize) '-zsc' rzsc '-window' int2str(window(2)) '-ctype-' ctype '.svg'],f2);
  close all
end

%  keyboard

