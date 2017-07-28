function RippleMod_PrePost_2(pval,varargin)


%RippleModPrePost - Plots ColorMaps for Ripple CCGS for a group of cells in pre and post-learning sleep epochs.
%
%  USAGE
%
%    RippleMod_PrePost_2(pval,<options>)
%
%    pval	    pvalue for ripple modulation
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%    structure          rat/sess/shank list for structure
%    CellList		CellList : rat/sess/shank
%    savefig		'on'/'off'
%    =========================================================================
%
%  OUTPUT
%
%   Figures
%
%  NOTE
%
%  SEE
%
%    See also
%
% Gabrielle Girardeau, June 2016
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
CellList=[];
savefig='off';

% Check number of inputs
if nargin < 1,
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
    case 'celllist',
      CellList = varargin{i+1};	
    case 'savefig'
      savefig = varargin{i+1};
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

load('/media/Data-01/All-Rats/SpikeParameters.mat');
load('/media/Data-01/All-Rats/PoissonRippleMod.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/AllRats-AllRippleCCGs.mat');

cellgroup = input('CellGroup?','s');

% Index the cells that are part of cell list
select.prepost=(ismember(idxall.prepost(:,1:4),CellList,'rows'));
select.all=(ismember(idxall.all(:,1:4),CellList,'rows'));

%Restrict poissonripplemod to all prepost index cells
poisson.prepost=poissonripplemod(ismember(poissonripplemod(:,1:4),idxall.prepost(:,1:4),'rows'),:);
poisson.all=poissonripplemod(ismember(poissonripplemod(:,1:4),idxall.all(:,1:4),'rows'),:);

ismod.prepost=poisson.prepost(:,6)<=pval/2|poisson.prepost(:,7)<=pval/2;
isInc.prepost=poisson.prepost(:,6)<=pval/2;
isDec.prepost=poisson.prepost(:,7)<=pval/2;

ismod.all=poisson.all(:,6)<=pval/2|poisson.all(:,7)<=pval/2;
isInc.all=poisson.all(:,6)<=pval/2;
isDec.all=poisson.all(:,7)<=pval/2;

%%% Select CCGs for Listed Cells & sig inc.
suballccgs.all.t1.inc=allccg.gain1(select.all&isInc.all,:);
suballccgs.pre.t1.inc=allccg.pre.gain1(select.prepost&isInc.prepost,:);
suballccgs.post.t1.inc=allccg.post.gain1(select.prepost&isInc.prepost,:);

suballccgs.all.t2.inc=allccg.gain2(select.all&isInc.all,:);
suballccgs.pre.t2.inc=allccg.pre.gain2(select.prepost&isInc.prepost,:);
suballccgs.post.t2.inc=allccg.post.gain2(select.prepost&isInc.prepost,:);


suballccgs.all.t1.dec=allccg.gain1(select.all&isDec.all,:);
suballccgs.pre.t1.dec=allccg.pre.gain1(select.prepost&isDec.prepost,:);
suballccgs.post.t1.dec=allccg.post.gain1(select.prepost&isDec.prepost,:);

suballccgs.all.t2.dec=allccg.gain2(select.all&isDec.all,:);
suballccgs.pre.t2.dec=allccg.pre.gain2(select.prepost&isDec.prepost,:);
suballccgs.post.t2.dec=allccg.post.gain2(select.prepost&isDec.prepost,:);

%% Eliminate cells with no data (for prepost comparison)
elim.t1=sum(suballccgs.pre.t1.inc,2)==0|sum(suballccgs.post.t1.inc,2)==0|isnan(suballccgs.pre.t1.inc(:,1))|isnan(suballccgs.post.t1.inc(:,1));
suballccgs.pre.t1.inc(elim.t1,:)=[];
suballccgs.post.t1.inc(elim.t1,:)=[];
suballccgs.pre.t2.inc(elim.t1,:)=[];
suballccgs.post.t2.inc(elim.t1,:)=[];

Delim.t1=sum(suballccgs.pre.t1.inc,2)==0|sum(suballccgs.post.t1.inc,2)==0|isnan(suballccgs.pre.t1.inc(:,1))|isnan(suballccgs.post.t1.inc(:,1));
suballccgs.pre.t1.dec(Delim.t1,:)=[];
suballccgs.post.t1.dec(Delim.t1,:)=[];
suballccgs.pre.t2.dec(Delim.t1,:)=[];
suballccgs.post.t2.dec(Delim.t1,:)=[];


%%% Stats
for i=1:size(suballccgs.pre.t2.inc,1)
  smoothccg.pre(i,:)=Smooth(suballccgs.pre.t2.inc(i,:),2);
  smoothccg.post(i,:)=Smooth(suballccgs.post.t2.inc(i,:),2);
%    meanpre(i)=mean(smoothccg.pre(i,190:211)); %mean 200ms window around ripple peak
%    meanpost(i)=mean(smoothccg.post(i,190:211));
  meanpre(i)=mean(smoothccg.pre(i,175:226)); %mean 500ms window around ripples peak
  meanpost(i)=mean(smoothccg.post(i,175:226));
end

%  for i=1:401
%    [ttp(i),h]=signrank(smoothccg.pre(:,i),smoothccg.post(:,i));
%  end
%  sig=ttp<0.01;
%  sigint=ToIntervals(allccg.t2,sig)
%  if size(sigint,2)==1;
%    sigint=sigint'
%  end


a=ones(length(meanpre),1);
figure;hold on;
plot([meanpre;meanpost]);
xlim([0 3])

[h3,p,stats]=signrank(meanpre,meanpost,'tail','left')
size(smoothccg.pre)

figure;
subplot(2,2,1)
idx=OrderedImagesc(suballccgs.pre.t2.inc,'newfig','off','x',allccg.t2,'norm','on');
subplot(2,2,2)
OrderedImagesc(suballccgs.post.t2.inc,'newfig','off','x',allccg.t2,'idx',idx,'norm','on');
subplot(2,2,3);hold on;
plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.inc),2),'b');
plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.inc),2)+Smooth(sem(suballccgs.pre.t2.inc),2),'b:');
plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.inc),2)-Smooth(sem(suballccgs.pre.t2.inc),2),'b:');
plot(allccg.t2,Smooth(mean(suballccgs.post.t2.inc),2),'r');
ylim([1 2]);
plot(allccg.t2,Smooth(mean(suballccgs.post.t2.inc),2)+Smooth(sem(suballccgs.post.t2.inc),2),'r:');
plot(allccg.t2,Smooth(mean(suballccgs.post.t2.inc),2)-Smooth(sem(suballccgs.post.t2.inc),2),'r:');

suptitle(['Ripple Mod + Cells - Norm and mean gains - Pre/Post - ' cellgroup])

subidx.inc=idxall.prepost(select.prepost&isInc.prepost,:);
subidx.inc(elim.t1,:)=[];
subidx.inc(idx,:);

%%%%%%%%%%%%%%%%% Decreasing cells
%  figure;
%  subplot(2,2,1)
%  idx=OrderedImagesc(suballccgs.pre.t2.dec,'newfig','off','x',allccg.t2,'norm','on');
%  subplot(2,2,2)
%  OrderedImagesc(suballccgs.post.t2.dec,'newfig','off','x',allccg.t2,'idx',idx,'norm','on');
%  subplot(2,2,3);hold on;
%  plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.dec),2),'b');
%  plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.dec),2)+Smooth(sem(suballccgs.pre.t2.dec),2),'b:');
%  plot(allccg.t2,Smooth(mean(suballccgs.pre.t2.dec),2)-Smooth(sem(suballccgs.pre.t2.dec),2),'b:');
%  plot(allccg.t2,Smooth(mean(suballccgs.post.t2.dec),2),'r');
%  ylim([0.5 1.1]);
%  plot(allccg.t2,Smooth(mean(suballccgs.post.t2.dec),2)+Smooth(sem(suballccgs.post.t2.dec),2),'r:');
%  plot(allccg.t2,Smooth(mean(suballccgs.post.t2.dec),2)-Smooth(sem(suballccgs.post.t2.dec),2),'r:');
%  suptitle(['Ripple Mod + Cells - Norm and mean gains - Pre/Post - ' cellgroup])
%  
%  subidx.dec=idxall.prepost(select.prepost&isDec.prepost,:);
%  subidx.dec(Delim.t1,:)=[];
%  subidx.dec(idx,:);


if strcmp(savefig,'on')
 cd ('/media/Data-01/All-Rats/AllRats-RippleMod-PrePost')
 plot2svg(['PrePostRippleGains-RippleModCells-CCGs-' cellgroup '.svg'],gcf);
end

