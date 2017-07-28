function [EV,REV,CorrMatrix] = ExplainedVariance_Ultimate(session,structure,pre,post,varargin)

% Explained Variance_Ultimate - Calculates the explained variance and reverse explained variances for pre/run/post sessions (Kudrimoti 1999, Lansink 2009)
%
%  USAGE
%
%    [EV,REV,CorrMatrix] = ExplainedVariance_Ultimate(session,structure,pre,post,varargin)
%
%    session		Complete path to session
%    structure		structure name (ex : 'BLA')
%    pre		subsession number for presleep
%    post		subsession number for postsleep
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'binsize'		size of bin for binning spike trains in seconds. Default : 0.05s
%     'savefig'		save figure (Default 'off')
%     'savevar'		save variable (Default 'on')		
%    =========================================================================
%
%  OUTPUT
%
%    EV           	Explained Variance
%    REV            	Reverse Explained variance
%    CorrMatrix		Correlation matrices for pre/run/post (CorrMatrix.pre .run .post)
%
%  NOTE
%
%
%  SEE
%
%    See also : binspikes, corrcoef, corr
%
% December 2015, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
binsize = 0.05;
savefig = 'off';
savevar = 'off';

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
		case 'brainstate',
			brainstate = lower(varargin{i+1});
		case 'binsize',
			binsize = (varargin{i+1});
		case 'savevar',
			savevar = (varargin{i+1});
		case 'savefig',
			savefig = (varargin{i+1});
		case 'ripplemod',
			ripplemod = (varargin{i+1});
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

load('/media/Data-01/All-Rats/sessionindexing.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
ratsess=ratsessionindex(strcmp([session '/'],xmlpath),1:2); %[rat session]

cd(session);
xml=session(end-13:end);
SetCurrentSession([session '/' xml]);

load([xml '-PoissonRippleMod.mat']);
load('States.mat');
load('runintervals.mat');
load(['/media/Data-01/All-Rats/Structures/' structure '.mat']);
load(['/media/Data-01/All-Rats/Structures/Hpc.mat']);
struc=eval(structure);

% Get pre/post sleep intervals
presleep=RunIntervals(pre);
postsleep=RunIntervals(post);
run=runintervals(2,:);

% Get SWS ripple intervals
ripples=GetRippleEvents;
swsripples=Restrict(ripples,sws);
swsrippleintervals=[swsripples(:,1) swsripples(:,2)];

swsrippleintervals_pre=Restrict(swsrippleintervals,presleep);
swsrippleintervals_post=Restrict(swsrippleintervals,postsleep);

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

% Ripple Modulation for structure cells Pval limit hardcoded = 0.001 for ripple mod.
STmod=poissonripplemod(ismember(poissonripplemod(:,1:2),STcells,'rows'),4:6); %local poisson variable (no rat/session columns)
is.up=STmod(:,1)<0.001/2;
is.down=STmod(:,2)<0.001/2;
is.mod=is.up|is.down;

endlimit=runintervals(end,end);
limits=[0 endlimit];

% Remove extraspikes for session with rec after run-post (ex : Rat08-20130710)
HPspikes(HPspikes(:,1)>=limits(2),:)=[];
STspikes(STspikes(:,1)>=limits(2),:)=[];

% Bin spikes
[HPtrains,bins]=SpikeTrain(HPspikes(:,[1 4]),binsize,limits);
[STtrains,bins]=SpikeTrain(STspikes(:,[1 4]),binsize,limits);

nbins20min=1200/binsize;

is.pre.sws=InIntervals(bins,presleep)&InIntervals(bins,sws);
is.post.sws=InIntervals(bins,postsleep)&InIntervals(bins,sws);
is.pre.rem=InIntervals(bins,presleep)&InIntervals(bins,Rem);
is.post.rem=InIntervals(bins,postsleep)&InIntervals(bins,Rem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate correlation matrices

%%%%%%%%%% SWS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
duration.pre.sws=sum(is.pre.sws)*binsize;
duration.post.sws=sum(is.post.sws)*binsize;
postswsHPtrains=HPtrains(is.post.sws,:);
postswsSTtrains=STtrains(is.post.sws,:);
  %%%%%%% All celltypes
    %%%% All modulation types
    [CorrMatrix.pre.sws.all.Rall, PvalMatrix.pre.sws.all.Rall]=corr(HPtrains(is.pre.sws,:),STtrains(is.pre.sws,:));
    [CorrMatrix.post.sws.all.Rall, PvalMatrix.post.sws.all.Rall]=corr(HPtrains(is.post.sws,:),STtrains(is.post.sws,:));
    [CorrMatrix.post_1.sws.all.Rall, PvalMatrix.post_1.sws.all.Rall]=corr(postswsHPtrains(1:nbins20min,:),postswsSTtrains(1:nbins20min,:));
    if size(postswsHPtrains,1)<2*nbins20min
      [CorrMatrix.post_2.sws.all.Rall, PvalMatrix.post_2.sws.all.Rall]=corr(postswsHPtrains(nbins20min+1:end,:),postswsSTtrains(nbins20min+1:end,:)); 
      CorrMatrix.post_3.sws.all.Rall=[];
      PvalMatrix.post_3.sws.all.Rall=[];
    else
      [CorrMatrix.post_2.sws.all.Rall, PvalMatrix.post_2.sws.all.Rall]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,:),postswsSTtrains(nbins20min+1:2*nbins20min,:));
      if size(postswsHPtrains,1)<3*nbins20min
	[CorrMatrix.post_3.sws.all.Rall, PvalMatrix.post_3.sws.all.Rall]=corr(postswsHPtrains(2*nbins20min+1:end,:),postswsSTtrains(2*nbins20min+1:end,:)); 
      else
	[CorrMatrix.post_3.sws.all.Rall, PvalMatrix.post_3.sws.all.Rall]=corr(postswsHPtrains(2*nbins20min+1:3*nbins20min,:),postswsSTtrains(2*nbins20min+1:3*nbins20min,:)); 
      end
    end

    %%%% Ripple Mod (all) vs non ripple mod
    [CorrMatrix.pre.sws.all.Rmod, PvalMatrix.pre.sws.all.Rmod]=corr(HPtrains(is.pre.sws,:),STtrains(is.pre.sws,is.mod));
    [CorrMatrix.post.sws.all.Rmod, PvalMatrix.post.sws.all.Rmod]=corr(HPtrains(is.post.sws,:),STtrains(is.post.sws,is.mod));
    [CorrMatrix.pre.sws.all.Rnomod, PvalMatrix.pre.sws.all.Rnomod]=corr(HPtrains(is.pre.sws,:),STtrains(is.pre.sws,~is.mod));
    [CorrMatrix.post.sws.all.Rnomod, PvalMatrix.post.sws.all.Rnomod]=corr(HPtrains(is.post.sws,:),STtrains(is.post.sws,~is.mod));    

    [CorrMatrix.post_1.sws.all.Rmod, PvalMatrix.post_1.sws.all.Rmod]=corr(postswsHPtrains(1:nbins20min,:),postswsSTtrains(1:nbins20min,is.mod));
    [CorrMatrix.post_1.sws.all.Rnomod, PvalMatrix.post_1.sws.all.Rnomod]=corr(postswsHPtrains(1:nbins20min,:),postswsSTtrains(1:nbins20min,~is.mod));
    if size(postswsHPtrains,1)<2*nbins20min
      [CorrMatrix.post_2.sws.all.Rmod, PvalMatrix.post_2.sws.all.Rmod]=corr(postswsHPtrains(nbins20min+1:end,:),postswsSTtrains(nbins20min+1:end,is.mod)); 
      CorrMatrix.post_3.sws.all.Rmod=[];
      PvalMatrix.post_3.sws.all.Rmod=[];
      [CorrMatrix.post_2.sws.all.Rnomod, PvalMatrix.post_2.sws.all.Rnomod]=corr(postswsHPtrains(nbins20min+1:end,:),postswsSTtrains(nbins20min+1:end,~is.mod)); 
      CorrMatrix.post_3.sws.all.Rnomod=[];
      PvalMatrix.post_3.sws.all.Rnomod=[];      
    else
      [CorrMatrix.post_2.sws.all.Rnomod, PvalMatrix.post_2.sws.all.Rnomod]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,:),postswsSTtrains(nbins20min+1:2*nbins20min,~is.mod));
      [CorrMatrix.post_2.sws.all.Rmod, PvalMatrix.post_2.sws.all.Rmod]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,:),postswsSTtrains(nbins20min+1:2*nbins20min,is.mod));
      if size(postswsHPtrains,1)<3*nbins20min
	[CorrMatrix.post_3.sws.all.Rmod, PvalMatrix.post_3.sws.all.Rmod]=corr(postswsHPtrains(2*nbins20min+1:end,:),postswsSTtrains(2*nbins20min+1:end,is.mod)); 
	[CorrMatrix.post_3.sws.all.Rnomod, PvalMatrix.post_3.sws.all.Rnomod]=corr(postswsHPtrains(2*nbins20min+1:end,:),postswsSTtrains(2*nbins20min+1:end,~is.mod)); 
      else
	[CorrMatrix.post_3.sws.all.Rmod, PvalMatrix.post_3.sws.all.Rmod]=corr(postswsHPtrains(2*nbins20min+1:3*nbins20min,:),postswsSTtrains(2*nbins20min+1:3*nbins20min,is.mod)); 
	[CorrMatrix.post_3.sws.all.Rnomod, PvalMatrix.post_3.sws.all.Rnomod]=corr(postswsHPtrains(2*nbins20min+1:end,:),postswsSTtrains(2*nbins20min+1:end,~is.mod)); 
      end
    end
    
  %%%%%%% Pyramidal cells only
    %%%% All modulation types
    [CorrMatrix.pre.sws.pyr.Rall, PvalMatrix.pre.sws.pyr.Rall]=corr(HPtrains(is.pre.sws,HPispyr),STtrains(is.pre.sws,STispyr));
    [CorrMatrix.post.sws.pyr.Rall, PvalMatrix.post.sws.pyr.Rall]=corr(HPtrains(is.post.sws,HPispyr),STtrains(is.post.sws,STispyr));

    [CorrMatrix.post_1.sws.pyr.Rall, PvalMatrix.post_1.sws.pyr.Rall]=corr(postswsHPtrains(1:nbins20min,HPispyr),postswsSTtrains(1:nbins20min,STispyr));    
    if size(postswsHPtrains,1)<2*nbins20min
      [CorrMatrix.post_2.sws.pyr.Rall, PvalMatrix.post_2.sws.pyr.Rall]=corr(postswsHPtrains(nbins20min+1:end,HPispyr),postswsSTtrains(nbins20min+1:end,STispyr)); 
      CorrMatrix.post_3.sws.pyr.Rall=[];
      PvalMatrix.post_3.sws.pyr.Rall=[];
    else
      [CorrMatrix.post_2.sws.pyr.Rall, PvalMatrix.post_2.sws.pyr.Rall]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,HPispyr),postswsSTtrains(nbins20min+1:2*nbins20min,STispyr));
      if size(postswsHPtrains,1)<3*nbins20min
	[CorrMatrix.post_3.sws.pyr.Rall, PvalMatrix.post_3.sws.pyr.Rall]=corr(postswsHPtrains(2*nbins20min+1:end,HPispyr),postswsSTtrains(2*nbins20min+1:end,STispyr)); 
      else
	[CorrMatrix.post_3.sws.pyr.Rall, PvalMatrix.post_3.sws.pyr.Rall]=corr(postswsHPtrains(2*nbins20min+1:3*nbins20min,HPispyr),postswsSTtrains(2*nbins20min+1:3*nbins20min,STispyr)); 
      end
    end

    %%%% Ripple Mod (all) vs non ripple mod
    [CorrMatrix.pre.sws.pyr.Rmod, PvalMatrix.pre.sws.pyr.Rmod]=corr(HPtrains(is.pre.sws,HPispyr),STtrains(is.pre.sws,is.mod&STispyr));
    [CorrMatrix.post.sws.pyr.Rmod, PvalMatrix.post.sws.pyr.Rmod]=corr(HPtrains(is.post.sws,HPispyr),STtrains(is.post.sws,is.mod&STispyr));
    [CorrMatrix.pre.sws.pyr.Rnomod, PvalMatrix.pre.sws.pyr.Rnomod]=corr(HPtrains(is.pre.sws,HPispyr),STtrains(is.pre.sws,~is.mod&STispyr));
    [CorrMatrix.post.sws.pyr.Rnomod, PvalMatrix.post.sws.pyr.Rnomod]=corr(HPtrains(is.post.sws,HPispyr),STtrains(is.post.sws,~is.mod&STispyr));    

    [CorrMatrix.post_1.sws.pyr.Rmod, PvalMatrix.post_1.sws.pyr.Rmod]=corr(postswsHPtrains(1:nbins20min,HPispyr),postswsSTtrains(1:nbins20min,is.mod&STispyr));
    [CorrMatrix.post_1.sws.pyr.Rnomod, PvalMatrix.post_1.sws.pyr.Rnomod]=corr(postswsHPtrains(1:nbins20min,HPispyr),postswsSTtrains(1:nbins20min,~is.mod&STispyr));
    if size(postswsHPtrains,1)<2*nbins20min
      [CorrMatrix.post_2.sws.pyr.Rnomod, PvalMatrix.post_2.sws.pyr.Rnomod]=corr(postswsHPtrains(nbins20min+1:end,HPispyr),postswsSTtrains(nbins20min+1:end,~is.mod&STispyr));
      [CorrMatrix.post_2.sws.pyr.Rmod, PvalMatrix.post_2.sws.pyr.Rmod]=corr(postswsHPtrains(nbins20min+1:end,HPispyr),postswsSTtrains(nbins20min+1:end,is.mod&STispyr));
      CorrMatrix.post_3.sws.pyr.Rmod=[];
      PvalMatrix.post_3.sws.pyr.Rmod=[];
      CorrMatrix.post_3.sws.pyr.Rnomod=[];
      PvalMatrix.post_3.sws.pyr.Rnomod=[];
    else 
      [CorrMatrix.post_2.sws.pyr.Rnomod, PvalMatrix.post_2.sws.pyr.Rnomod]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,HPispyr),postswsSTtrains(nbins20min+1:2*nbins20min,~is.mod&STispyr));
      [CorrMatrix.post_2.sws.pyr.Rmod, PvalMatrix.post_2.sws.pyr.Rmod]=corr(postswsHPtrains(nbins20min+1:2*nbins20min,HPispyr),postswsSTtrains(nbins20min+1:2*nbins20min,is.mod&STispyr));
      if size(postswsHPtrains,1)<3*nbins20min
	[CorrMatrix.post_3.sws.pyr.Rmod, PvalMatrix.post_3.sws.pyr.Rmod]=corr(postswsHPtrains(2*nbins20min+1:end,HPispyr),postswsSTtrains(2*nbins20min+1:end,is.mod&STispyr)); 
	[CorrMatrix.post_3.sws.pyr.Rnomod, PvalMatrix.post_3.sws.pyr.Rnomod]=corr(postswsHPtrains(2*nbins20min+1:end,HPispyr),postswsSTtrains(2*nbins20min+1:end,~is.mod&STispyr)); 
      else
	[CorrMatrix.post_3.sws.pyr.Rmod, PvalMatrix.post_3.sws.pyr.Rmod]=corr(postswsHPtrains(2*nbins20min+1:3*nbins20min,HPispyr),postswsSTtrains(2*nbins20min+1:3*nbins20min,is.mod&STispyr)); 
	[CorrMatrix.post_3.sws.pyr.Rnomod, PvalMatrix.post_3.sws.pyr.Rnomod]=corr(postswsHPtrains(2*nbins20min+1:end,HPispyr),postswsSTtrains(2*nbins20min+1:end,~is.mod&STispyr)); 
      end
    end
    
%%%%%%%%%%% REM %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration.pre.rem=sum(InIntervals(bins,presleep)&InIntervals(bins,Rem))*binsize;
duration.post.rem=sum(InIntervals(bins,postsleep)&InIntervals(bins,Rem))*binsize;
  %%%%%%%% AllCelltypes
    %%%%%%%%%% All ripple mod types
    [CorrMatrix.pre.rem.all.Rall, PvalMatrix.pre.rem.all.Rall]=corr(HPtrains(is.pre.rem,:),STtrains(is.pre.rem,:));
    [CorrMatrix.post.rem.all.Rall, PvalMatrix.post.rem.all.Rall]=corr(HPtrains(is.post.rem,:),STtrains(is.post.rem,:));
    
    %%%%%%%%%% Mod vs No mod
    [CorrMatrix.pre.rem.all.Rmod, PvalMatrix.pre.rem.all.Rmod]=corr(HPtrains(is.pre.rem,:),STtrains(is.pre.rem,is.mod));
    [CorrMatrix.post.rem.all.Rmod, PvalMatrix.post.rem.all.Rmod]=corr(HPtrains(is.post.rem,:),STtrains(is.post.rem,is.mod));
    [CorrMatrix.pre.rem.all.Rnomod, PvalMatrix.pre.rem.all.Rnomod]=corr(HPtrains(is.pre.rem,:),STtrains(is.pre.rem,~is.mod));
    [CorrMatrix.post.rem.all.Rnomod, PvalMatrix.post.rem.all.Rnomod]=corr(HPtrains(is.post.rem,:),STtrains(is.post.rem,~is.mod));
       
  %%%%%% Pyr Only
    %%%%%%%%% All Ripple mod types
    [CorrMatrix.pre.rem.pyr.Rall, PvalMatrix.pre.rem.pyr.Rall]=corr(HPtrains(is.pre.rem,HPispyr),STtrains(is.pre.rem,STispyr));
    [CorrMatrix.post.rem.pyr.Rall, PvalMatrix.post.rem.pyr.Rall]=corr(HPtrains(is.post.rem,HPispyr),STtrains(is.post.rem,STispyr));
    
    %%%%%%%%%% Mod vs No mod
    [CorrMatrix.pre.rem.pyr.Rmod, PvalMatrix.pre.rem.pyr.Rmod]=corr(HPtrains(is.pre.rem,HPispyr),STtrains(is.pre.rem,is.mod&STispyr));
    [CorrMatrix.post.rem.pyr.Rmod, PvalMatrix.post.rem.pyr.Rmod]=corr(HPtrains(is.post.rem,HPispyr),STtrains(is.post.rem,is.mod&STispyr));
    [CorrMatrix.pre.rem.pyr.Rnomod, PvalMatrix.pre.rem.pyr.Rnomod]=corr(HPtrains(is.pre.rem,HPispyr),STtrains(is.pre.rem,~is.mod&STispyr));
    [CorrMatrix.post.rem.pyr.Rnomod, PvalMatrix.post.rem.pyr.Rnomod]=corr(HPtrains(is.post.rem,HPispyr),STtrains(is.post.rem,~is.mod&STispyr));
       
%%%%%%%%%%%% Ripple In/Out %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is.inrip.pre=InIntervals(bins,[swsrippleintervals_pre(:,1)-0.05 swsrippleintervals_pre(:,2)]);
is.outrip.pre=~InIntervals(bins,[swsrippleintervals_pre(:,1)-0.05 swsrippleintervals_pre(:,2)]);
is.inrip.post=InIntervals(bins,[swsrippleintervals_post(:,1)-0.05 swsrippleintervals_post(:,2)]);
is.outrip.post=~InIntervals(bins,[swsrippleintervals_post(:,1)-0.05 swsrippleintervals_post(:,2)]);
duration.pre.rip=sum(InIntervals(bins,[swsrippleintervals_pre(:,1)-0.05 swsrippleintervals_pre(:,2)]))*binsize;
duration.post.rip=sum(InIntervals(bins,[swsrippleintervals_post(:,1)-0.05 swsrippleintervals_post(:,2)]))*binsize;
  %%%%% AllCellTypes
    %% All ripple modulation types
    [CorrMatrix.pre.rip_in.all.Rall ,PvalMatrix.pre.rip_in.all.Rall]=corr(HPtrains(is.inrip.pre,:),STtrains(is.inrip.pre,:));
    [CorrMatrix.post.rip_in.all.Rall, PvalMatrix.post.rip_in.all.Rall]=corr(HPtrains(is.inrip.post,:),STtrains(is.inrip.post,:));
    [CorrMatrix.pre.rip_out.all.Rall, PvalMatrix.pre.rip_out.all.Rall]=corr(HPtrains(~is.inrip.pre&is.pre.sws,:),STtrains(~is.inrip.pre&is.pre.sws,:));
    [CorrMatrix.post.rip_out.all.Rall, PvalMatrix.post.rip_out.all.Rall]=corr(HPtrains(~is.inrip.post&is.post.sws,:),STtrains(~is.inrip.post&is.post.sws,:));
    %% Mod vs No mod
    [CorrMatrix.pre.rip_in.all.Rmod ,PvalMatrix.pre.rip_in.all.Rmod]=corr(HPtrains(is.inrip.pre,:),STtrains(is.inrip.pre,is.mod));
    [CorrMatrix.post.rip_in.all.Rmod, PvalMatrix.post.rip_in.all.Rmod]=corr(HPtrains(is.inrip.post,:),STtrains(is.inrip.post,is.mod));
    [CorrMatrix.pre.rip_out.all.Rmod, PvalMatrix.pre.rip_out.all.Rmod]=corr(HPtrains(~is.inrip.pre&is.pre.sws,:),STtrains(~is.inrip.pre&is.pre.sws,is.mod));
    [CorrMatrix.post.rip_out.all.Rmod, PvalMatrix.post.rip_out.all.Rmod]=corr(HPtrains(~is.inrip.post&is.post.sws,:),STtrains(~is.inrip.post&is.post.sws,is.mod));

    [CorrMatrix.pre.rip_in.all.Rnomod ,PvalMatrix.pre.rip_in.all.Rnomod]=corr(HPtrains(is.inrip.pre,:),STtrains(is.inrip.pre,~is.mod));
    [CorrMatrix.post.rip_in.all.Rnomod, PvalMatrix.post.rip_in.all.Rnomod]=corr(HPtrains(is.inrip.post,:),STtrains(is.inrip.post,~is.mod));
    [CorrMatrix.pre.rip_out.all.Rnomod, PvalMatrix.pre.rip_out.all.Rnomod]=corr(HPtrains(~is.inrip.pre&is.pre.sws,:),STtrains(~is.inrip.pre&is.pre.sws,~is.mod));
    [CorrMatrix.post.rip_out.all.Rnomod, PvalMatrix.post.rip_out.all.Rmod]=corr(HPtrains(~is.inrip.post&is.post.sws,:),STtrains(~is.inrip.post&is.post.sws,~is.mod));

  %%%%%%%% Pyr only
    %% All ripple modulation types
    [CorrMatrix.pre.rip_in.pyr.Rall ,PvalMatrix.pre.rip_in.pyr.Rall]=corr(HPtrains(is.inrip.pre,HPispyr),STtrains(is.inrip.pre,STispyr));
    [CorrMatrix.post.rip_in.pyr.Rall, PvalMatrix.post.rip_in.pyr.Rall]=corr(HPtrains(is.inrip.post,HPispyr),STtrains(is.inrip.post,STispyr));
    [CorrMatrix.pre.rip_out.pyr.Rall, PvalMatrix.pre.rip_out.pyr.Rall]=corr(HPtrains(~is.inrip.pre&is.pre.sws,HPispyr),STtrains(~is.inrip.pre&is.pre.sws,STispyr));
    [CorrMatrix.post.rip_out.pyr.Rall, PvalMatrix.post.rip_out.pyr.Rall]=corr(HPtrains(~is.inrip.post&is.post.sws,HPispyr),STtrains(~is.inrip.post&is.post.sws,STispyr));
    %% Mod vs No mod
    [CorrMatrix.pre.rip_in.pyr.Rmod ,PvalMatrix.pre.rip_in.pyr.Rmod]=corr(HPtrains(is.inrip.pre,HPispyr),STtrains(is.inrip.pre,is.mod&STispyr));
    [CorrMatrix.post.rip_in.pyr.Rmod, PvalMatrix.post.rip_in.pyr.Rmod]=corr(HPtrains(is.inrip.post,HPispyr),STtrains(is.inrip.post,is.mod&STispyr));
    [CorrMatrix.pre.rip_out.pyr.Rmod, PvalMatrix.pre.rip_out.pyr.Rmod]=corr(HPtrains(~is.inrip.pre&is.pre.sws,HPispyr),STtrains(~is.inrip.pre&is.pre.sws,is.mod&STispyr));
    [CorrMatrix.post.rip_out.pyr.Rmod, PvalMatrix.post.rip_out.pyr.Rmod]=corr(HPtrains(~is.inrip.post&is.post.sws,HPispyr),STtrains(~is.inrip.post&is.post.sws,is.mod&STispyr));

    [CorrMatrix.pre.rip_in.pyr.Rnomod ,PvalMatrix.pre.rip_in.pyr.Rnomod]=corr(HPtrains(is.inrip.pre,HPispyr),STtrains(is.inrip.pre,~is.mod&STispyr));
    [CorrMatrix.post.rip_in.pyr.Rnomod, PvalMatrix.post.rip_in.pyr.Rnomod]=corr(HPtrains(is.inrip.post,HPispyr),STtrains(is.inrip.post,~is.mod&STispyr));
    [CorrMatrix.pre.rip_out.pyr.Rnomod, PvalMatrix.pre.rip_out.pyr.Rnomod]=corr(HPtrains(~is.inrip.pre&is.pre.sws,HPispyr),STtrains(~is.inrip.pre&is.pre.sws,~is.mod&STispyr));
    [CorrMatrix.post.rip_out.pyr.Rnomod, PvalMatrix.post.rip_out.pyr.Rmod]=corr(HPtrains(~is.inrip.post&is.post.sws,HPispyr),STtrains(~is.inrip.post&is.post.sws,~is.mod&STispyr));

%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is.run=InIntervals(bins,runintervals(2,:));
  %%%%%% Allcell type
    % All Ripple mod
    [CorrMatrix.run.all.Rall, PvalMatrix.run.all.Rall]=corr(HPtrains(is.run,:),STtrains(is.run,:));   
    % Ripple mod vs Non ripplemod
    [CorrMatrix.run.all.Rmod, PvalMatrix.run.all.Rmod]=corr(HPtrains(is.run,:),STtrains(is.run,is.mod));   
    [CorrMatrix.run.all.Rnomod, PvalMatrix.run.all.Rnomod]=corr(HPtrains(is.run,:),STtrains(is.run,~is.mod)); 

  %%%%%% Pyramidals Only
    % All Ripple mod
    [CorrMatrix.run.pyr.Rall, PvalMatrix.run.pyr.Rall]=corr(HPtrains(is.run,HPispyr),STtrains(is.run,STispyr));   
    % Ripple mod vs Non ripplemod
    [CorrMatrix.run.pyr.Rmod, PvalMatrix.run.pyr.Rmod]=corr(HPtrains(is.run,HPispyr),STtrains(is.run,is.mod&STispyr));   
    [CorrMatrix.run.pyr.Rnomod, PvalMatrix.run.pyr.Rnomod]=corr(HPtrains(is.run,HPispyr),STtrains(is.run,~is.mod&STispyr)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation Coeffs.

%%%%%%%%%%%%% REM %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% All cell types
    %%%% All Ripple modulation types
    r_runpost.rem.all.Rall=corrcoef(CorrMatrix.post.rem.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost.rem.all.Rall=r_runpost.rem.all.Rall(1,2);
    r_runpre.rem.all.Rall=corrcoef(CorrMatrix.pre.rem.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpre.rem.all.Rall=r_runpre.rem.all.Rall(1,2);
    r_prepost.rem.all.Rall=corrcoef(CorrMatrix.pre.rem.all.Rall,CorrMatrix.post.rem.all.Rall,'rows','complete');r_prepost.rem.all.Rall=r_prepost.rem.all.Rall(1,2);
    EV.rem.all.Rall=((r_runpost.rem.all.Rall-r_runpre.rem.all.Rall*r_prepost.rem.all.Rall)/sqrt((1-r_runpre.rem.all.Rall^2)*(1-r_prepost.rem.all.Rall^2)))^2;
    REV.rem.all.Rall=((r_runpre.rem.all.Rall-r_runpost.rem.all.Rall*r_prepost.rem.all.Rall)/sqrt((1-r_runpost.rem.all.Rall^2)*(1-r_prepost.rem.all.Rall^2)))^2;
    
    %%%% Ripple mod vs nonripple mod
      %%% Ripple mod all
      r_runpost.rem.all.Rmod=corrcoef(CorrMatrix.post.rem.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost.rem.all.Rmod=r_runpost.rem.all.Rmod(1,2);
      r_runpre.rem.all.Rmod=corrcoef(CorrMatrix.pre.rem.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpre.rem.all.Rmod=r_runpre.rem.all.Rmod(1,2);
      r_prepost.rem.all.Rmod=corrcoef(CorrMatrix.pre.rem.all.Rmod,CorrMatrix.post.rem.all.Rmod,'rows','complete');r_prepost.rem.all.Rmod=r_prepost.rem.all.Rmod(1,2);
      EV.rem.all.Rmod=((r_runpost.rem.all.Rmod-r_runpre.rem.all.Rmod*r_prepost.rem.all.Rmod)/sqrt((1-r_runpre.rem.all.Rmod^2)*(1-r_prepost.rem.all.Rmod^2)))^2;
      REV.rem.all.Rmod=((r_runpre.rem.all.Rmod-r_runpost.rem.all.Rmod*r_prepost.rem.all.Rmod)/sqrt((1-r_runpost.rem.all.Rmod^2)*(1-r_prepost.rem.all.Rmod^2)))^2;
      
      %%% Ripple No mod
      r_runpost.rem.all.Rnomod=corrcoef(CorrMatrix.post.rem.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost.rem.all.Rnomod=r_runpost.rem.all.Rnomod(1,2);
      r_runpre.rem.all.Rnomod=corrcoef(CorrMatrix.pre.rem.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpre.rem.all.Rnomod=r_runpre.rem.all.Rnomod(1,2);
      r_prepost.rem.all.Rnomod=corrcoef(CorrMatrix.pre.rem.all.Rnomod,CorrMatrix.post.rem.all.Rnomod,'rows','complete');r_prepost.rem.all.Rnomod=r_prepost.rem.all.Rnomod(1,2);
      EV.rem.all.Rnomod=((r_runpost.rem.all.Rnomod-r_runpre.rem.all.Rnomod*r_prepost.rem.all.Rnomod)/sqrt((1-r_runpre.rem.all.Rnomod^2)*(1-r_prepost.rem.all.Rnomod^2)))^2;
      REV.rem.all.Rnomod=((r_runpre.rem.all.Rnomod-r_runpost.rem.all.Rnomod*r_prepost.rem.all.Rnomod)/sqrt((1-r_runpost.rem.all.Rnomod^2)*(1-r_prepost.rem.all.Rnomod^2)))^2;
          
  %%%% Pyr Only
    %%%%% All ripple mod types
    r_runpost.rem.pyr.Rall=corrcoef(CorrMatrix.post.rem.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost.rem.pyr.Rall=r_runpost.rem.pyr.Rall(1,2);
    r_runpre.rem.pyr.Rall=corrcoef(CorrMatrix.pre.rem.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpre.rem.pyr.Rall=r_runpre.rem.pyr.Rall(1,2);
    r_prepost.rem.pyr.Rall=corrcoef(CorrMatrix.pre.rem.pyr.Rall,CorrMatrix.post.rem.pyr.Rall,'rows','complete');r_prepost.rem.pyr.Rall=r_prepost.rem.pyr.Rall(1,2);
    EV.rem.pyr.Rall=((r_runpost.rem.pyr.Rall-r_runpre.rem.pyr.Rall*r_prepost.rem.pyr.Rall)/sqrt((1-r_runpre.rem.pyr.Rall^2)*(1-r_prepost.rem.pyr.Rall^2)))^2;
    REV.rem.pyr.Rall=((r_runpre.rem.pyr.Rall-r_runpost.rem.pyr.Rall*r_prepost.rem.pyr.Rall)/sqrt((1-r_runpost.rem.pyr.Rall^2)*(1-r_prepost.rem.pyr.Rall^2)))^2;

    %%%% Ripple mod vs nonripple mod
      %%% Ripple mod all
      r_runpost.rem.pyr.Rmod=corrcoef(CorrMatrix.post.rem.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost.rem.pyr.Rmod=r_runpost.rem.pyr.Rmod(1,2);
      r_runpre.rem.pyr.Rmod=corrcoef(CorrMatrix.pre.rem.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpre.rem.pyr.Rmod=r_runpre.rem.pyr.Rmod(1,2);
      r_prepost.rem.pyr.Rmod=corrcoef(CorrMatrix.pre.rem.pyr.Rmod,CorrMatrix.post.rem.pyr.Rmod,'rows','complete');r_prepost.rem.pyr.Rmod=r_prepost.rem.pyr.Rmod(1,2);
      EV.rem.pyr.Rmod=((r_runpost.rem.pyr.Rmod-r_runpre.rem.pyr.Rmod*r_prepost.rem.pyr.Rmod)/sqrt((1-r_runpre.rem.pyr.Rmod^2)*(1-r_prepost.rem.pyr.Rmod^2)))^2;
      REV.rem.pyr.Rmod=((r_runpre.rem.pyr.Rmod-r_runpost.rem.pyr.Rmod*r_prepost.rem.pyr.Rmod)/sqrt((1-r_runpost.rem.pyr.Rmod^2)*(1-r_prepost.rem.pyr.Rmod^2)))^2;
      
      %%% Ripple No mod
      r_runpost.rem.pyr.Rnomod=corrcoef(CorrMatrix.post.rem.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost.rem.pyr.Rnomod=r_runpost.rem.pyr.Rnomod(1,2);
      r_runpre.rem.pyr.Rnomod=corrcoef(CorrMatrix.pre.rem.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpre.rem.pyr.Rnomod=r_runpre.rem.pyr.Rnomod(1,2);
      r_prepost.rem.pyr.Rnomod=corrcoef(CorrMatrix.pre.rem.pyr.Rnomod,CorrMatrix.post.rem.pyr.Rnomod,'rows','complete');r_prepost.rem.pyr.Rnomod=r_prepost.rem.pyr.Rnomod(1,2);
      EV.rem.pyr.Rnomod=((r_runpost.rem.pyr.Rnomod-r_runpre.rem.pyr.Rnomod*r_prepost.rem.pyr.Rnomod)/sqrt((1-r_runpre.rem.pyr.Rnomod^2)*(1-r_prepost.rem.pyr.Rnomod^2)))^2;
      REV.rem.pyr.Rnomod=((r_runpre.rem.pyr.Rnomod-r_runpost.rem.pyr.Rnomod*r_prepost.rem.pyr.Rnomod)/sqrt((1-r_runpost.rem.pyr.Rnomod^2)*(1-r_prepost.rem.pyr.Rnomod^2)))^2;
              
%%%%%%%%%%%% SWS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% All cell types
    %%%%%% All ripple mod types
    r_runpost.sws.all.Rall=corrcoef(CorrMatrix.post.sws.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost.sws.all.Rall=r_runpost.sws.all.Rall(1,2);
    r_runpre.sws.all.Rall=corrcoef(CorrMatrix.pre.sws.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpre.sws.all.Rall=r_runpre.sws.all.Rall(1,2);
    r_prepost.sws.all.Rall=corrcoef(CorrMatrix.pre.sws.all.Rall,CorrMatrix.post.sws.all.Rall,'rows','complete');r_prepost.sws.all.Rall=r_prepost.sws.all.Rall(1,2);
    EV.sws.all.Rall=((r_runpost.sws.all.Rall-r_runpre.sws.all.Rall*r_prepost.sws.all.Rall)/sqrt((1-r_runpre.sws.all.Rall^2)*(1-r_prepost.sws.all.Rall^2)))^2;
    REV.sws.all.Rall=((r_runpre.sws.all.Rall-r_runpost.sws.all.Rall*r_prepost.sws.all.Rall)/sqrt((1-r_runpost.sws.all.Rall^2)*(1-r_prepost.sws.all.Rall^2)))^2;

    r_runpost_1.sws.all.Rall=corrcoef(CorrMatrix.post_1.sws.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost_1.sws.all.Rall=r_runpost_1.sws.all.Rall(1,2);
    r_runpost_2.sws.all.Rall=corrcoef(CorrMatrix.post_2.sws.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost_2.sws.all.Rall=r_runpost_2.sws.all.Rall(1,2);
    r_prepost_1.sws.all.Rall=corrcoef(CorrMatrix.pre.sws.all.Rall,CorrMatrix.post_1.sws.all.Rall,'rows','complete');r_prepost_1.sws.all.Rall=r_prepost_1.sws.all.Rall(1,2);
    r_prepost_2.sws.all.Rall=corrcoef(CorrMatrix.pre.sws.all.Rall,CorrMatrix.post_2.sws.all.Rall,'rows','complete');r_prepost_2.sws.all.Rall=r_prepost_2.sws.all.Rall(1,2);    
    EV.post1.sws.all.Rall=((r_runpost_1.sws.all.Rall-r_runpre.sws.all.Rall*r_prepost_1.sws.all.Rall)/sqrt((1-r_runpre.sws.all.Rall^2)*(1-r_prepost_1.sws.all.Rall^2)))^2;
    REV.post1.sws.all.Rall=((r_runpre.sws.all.Rall-r_runpost_1.sws.all.Rall*r_prepost_1.sws.all.Rall)/sqrt((1-r_runpost_1.sws.all.Rall^2)*(1-r_prepost_1.sws.all.Rall^2)))^2;
    EV.post2.sws.all.Rall=((r_runpost_2.sws.all.Rall-r_runpre.sws.all.Rall*r_prepost_2.sws.all.Rall)/sqrt((1-r_runpre.sws.all.Rall^2)*(1-r_prepost_2.sws.all.Rall^2)))^2;
    REV.post2.sws.all.Rall=((r_runpre.sws.all.Rall-r_runpost_2.sws.all.Rall*r_prepost_2.sws.all.Rall)/sqrt((1-r_runpost_2.sws.all.Rall^2)*(1-r_prepost_2.sws.all.Rall^2)))^2;  
    if ~isempty(CorrMatrix.post_3.sws.all.Rall)
      r_runpost_3.sws.all.Rall=corrcoef(CorrMatrix.post_3.sws.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost_3.sws.all.Rall=r_runpost_3.sws.all.Rall(1,2);
      r_prepost_3.sws.all.Rall=corrcoef(CorrMatrix.pre.sws.all.Rall,CorrMatrix.post_3.sws.all.Rall,'rows','complete');r_prepost_3.sws.all.Rall=r_prepost_3.sws.all.Rall(1,2);
      EV.post3.sws.all.Rall=((r_runpost_3.sws.all.Rall-r_runpre.sws.all.Rall*r_prepost_3.sws.all.Rall)/sqrt((1-r_runpre.sws.all.Rall^2)*(1-r_prepost_3.sws.all.Rall^2)))^2;
      REV.post3.sws.all.Rall=((r_runpre.sws.all.Rall-r_runpost_3.sws.all.Rall*r_prepost_3.sws.all.Rall)/sqrt((1-r_runpost_3.sws.all.Rall^2)*(1-r_prepost_3.sws.all.Rall^2)))^2;
    else
      r_runpost_3.sws.all.Rall=NaN;
      r_runpost_3.sws.all.Rall=NaN;
      EV.post3.sws.all.Rall=NaN;
      REV.post3.sws.all.Rall=NaN;
    end
    
    
    %%%% Ripple mod vs nonripple mod
      %%% Ripple mod all
      r_runpost.sws.all.Rmod=corrcoef(CorrMatrix.post.sws.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost.sws.all.Rmod=r_runpost.sws.all.Rmod(1,2);
      r_runpre.sws.all.Rmod=corrcoef(CorrMatrix.pre.sws.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpre.sws.all.Rmod=r_runpre.sws.all.Rmod(1,2);
      r_prepost.sws.all.Rmod=corrcoef(CorrMatrix.pre.sws.all.Rmod,CorrMatrix.post.sws.all.Rmod,'rows','complete');r_prepost.sws.all.Rmod=r_prepost.sws.all.Rmod(1,2);
      EV.sws.all.Rmod=((r_runpost.sws.all.Rmod-r_runpre.sws.all.Rmod*r_prepost.sws.all.Rmod)/sqrt((1-r_runpre.sws.all.Rmod^2)*(1-r_prepost.sws.all.Rmod^2)))^2;
      REV.sws.all.Rmod=((r_runpre.sws.all.Rmod-r_runpost.sws.all.Rmod*r_prepost.sws.all.Rmod)/sqrt((1-r_runpost.sws.all.Rmod^2)*(1-r_prepost.sws.all.Rmod^2)))^2;

      r_runpost_1.sws.all.Rmod=corrcoef(CorrMatrix.post_1.sws.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost_1.sws.all.Rmod=r_runpost_1.sws.all.Rmod(1,2);
      r_runpost_2.sws.all.Rmod=corrcoef(CorrMatrix.post_2.sws.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost_2.sws.all.Rmod=r_runpost_2.sws.all.Rmod(1,2);
      
      r_prepost_1.sws.all.Rmod=corrcoef(CorrMatrix.pre.sws.all.Rmod,CorrMatrix.post_1.sws.all.Rmod,'rows','complete');r_prepost_1.sws.all.Rmod=r_prepost_1.sws.all.Rmod(1,2);
      r_prepost_2.sws.all.Rmod=corrcoef(CorrMatrix.pre.sws.all.Rmod,CorrMatrix.post_2.sws.all.Rmod,'rows','complete');r_prepost_2.sws.all.Rmod=r_prepost_2.sws.all.Rmod(1,2);
      
      EV.post1.sws.all.Rmod=((r_runpost_1.sws.all.Rmod-r_runpre.sws.all.Rmod*r_prepost_1.sws.all.Rmod)/sqrt((1-r_runpre.sws.all.Rmod^2)*(1-r_prepost_1.sws.all.Rmod^2)))^2;
      REV.post1.sws.all.Rmod=((r_runpre.sws.all.Rmod-r_runpost_1.sws.all.Rmod*r_prepost_1.sws.all.Rmod)/sqrt((1-r_runpost_1.sws.all.Rmod^2)*(1-r_prepost_1.sws.all.Rmod^2)))^2;
      EV.post2.sws.all.Rmod=((r_runpost_2.sws.all.Rmod-r_runpre.sws.all.Rmod*r_prepost_2.sws.all.Rmod)/sqrt((1-r_runpre.sws.all.Rmod^2)*(1-r_prepost_2.sws.all.Rmod^2)))^2;
      REV.post2.sws.all.Rmod=((r_runpre.sws.all.Rmod-r_runpost_2.sws.all.Rmod*r_prepost_2.sws.all.Rmod)/sqrt((1-r_runpost_2.sws.all.Rmod^2)*(1-r_prepost_2.sws.all.Rmod^2)))^2;
      
      %%% Ripple nomod all
      r_runpost.sws.all.Rnomod=corrcoef(CorrMatrix.post.sws.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost.sws.all.Rnomod=r_runpost.sws.all.Rnomod(1,2);
      r_runpre.sws.all.Rnomod=corrcoef(CorrMatrix.pre.sws.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpre.sws.all.Rnomod=r_runpre.sws.all.Rnomod(1,2);
      r_prepost.sws.all.Rnomod=corrcoef(CorrMatrix.pre.sws.all.Rnomod,CorrMatrix.post.sws.all.Rnomod,'rows','complete');r_prepost.sws.all.Rnomod=r_prepost.sws.all.Rnomod(1,2);
      EV.sws.all.Rnomod=((r_runpost.sws.all.Rnomod-r_runpre.sws.all.Rnomod*r_prepost.sws.all.Rnomod)/sqrt((1-r_runpre.sws.all.Rnomod^2)*(1-r_prepost.sws.all.Rnomod^2)))^2;
      REV.sws.all.Rnomod=((r_runpre.sws.all.Rnomod-r_runpost.sws.all.Rnomod*r_prepost.sws.all.Rnomod)/sqrt((1-r_runpost.sws.all.Rnomod^2)*(1-r_prepost.sws.all.Rnomod^2)))^2;
      
      r_runpost_1.sws.all.Rnomod=corrcoef(CorrMatrix.post_1.sws.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost_1.sws.all.Rnomod=r_runpost_1.sws.all.Rnomod(1,2);
      r_runpost_2.sws.all.Rnomod=corrcoef(CorrMatrix.post_2.sws.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost_2.sws.all.Rnomod=r_runpost_2.sws.all.Rnomod(1,2);      
      r_prepost_1.sws.all.Rnomod=corrcoef(CorrMatrix.pre.sws.all.Rnomod,CorrMatrix.post_1.sws.all.Rnomod,'rows','complete');r_prepost_1.sws.all.Rnomod=r_prepost_1.sws.all.Rnomod(1,2);
      r_prepost_2.sws.all.Rnomod=corrcoef(CorrMatrix.pre.sws.all.Rnomod,CorrMatrix.post_2.sws.all.Rnomod,'rows','complete');r_prepost_2.sws.all.Rnomod=r_prepost_2.sws.all.Rnomod(1,2);    
      EV.post1.sws.all.Rnomod=((r_runpost_1.sws.all.Rnomod-r_runpre.sws.all.Rnomod*r_prepost_1.sws.all.Rnomod)/sqrt((1-r_runpre.sws.all.Rnomod^2)*(1-r_prepost_1.sws.all.Rnomod^2)))^2;
      REV.post1.sws.all.Rnomod=((r_runpre.sws.all.Rnomod-r_runpost_1.sws.all.Rnomod*r_prepost_1.sws.all.Rnomod)/sqrt((1-r_runpost_1.sws.all.Rnomod^2)*(1-r_prepost_1.sws.all.Rnomod^2)))^2;
      EV.post2.sws.all.Rnomod=((r_runpost_2.sws.all.Rnomod-r_runpre.sws.all.Rnomod*r_prepost_2.sws.all.Rnomod)/sqrt((1-r_runpre.sws.all.Rnomod^2)*(1-r_prepost_2.sws.all.Rnomod^2)))^2;
      REV.post2.sws.all.Rnomod=((r_runpre.sws.all.Rnomod-r_runpost_2.sws.all.Rnomod*r_prepost_2.sws.all.Rnomod)/sqrt((1-r_runpost_2.sws.all.Rnomod^2)*(1-r_prepost_2.sws.all.Rnomod^2)))^2     
      
      if ~isempty(CorrMatrix.post_3.sws.all.Rmod)
	r_runpost_3.sws.all.Rmod=corrcoef(CorrMatrix.post_3.sws.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost_3.sws.all.Rmod=r_runpost_3.sws.all.Rmod(1,2);
	r_prepost_3.sws.all.Rmod=corrcoef(CorrMatrix.pre.sws.all.Rmod,CorrMatrix.post_3.sws.all.Rmod,'rows','complete');r_prepost_3.sws.all.Rmod=r_prepost_3.sws.all.Rmod(1,2);
	EV.post3.sws.all.Rmod=((r_runpost_3.sws.all.Rmod-r_runpre.sws.all.Rmod*r_prepost_3.sws.all.Rmod)/sqrt((1-r_runpre.sws.all.Rmod^2)*(1-r_prepost_3.sws.all.Rmod^2)))^2;
	REV.post3.sws.all.Rmod=((r_runpre.sws.all.Rmod-r_runpost_3.sws.all.Rmod*r_prepost_3.sws.all.Rmod)/sqrt((1-r_runpost_3.sws.all.Rmod^2)*(1-r_prepost_3.sws.all.Rmod^2)))^2;
	r_runpost_3.sws.all.Rnomod=corrcoef(CorrMatrix.post_3.sws.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost_3.sws.all.Rnomod=r_runpost_3.sws.all.Rnomod(1,2);
	r_prepost_3.sws.all.Rnomod=corrcoef(CorrMatrix.pre.sws.all.Rnomod,CorrMatrix.post_3.sws.all.Rnomod,'rows','complete');r_prepost_3.sws.all.Rnomod=r_prepost_3.sws.all.Rnomod(1,2);
	EV.post3.sws.all.Rnomod=((r_runpost_3.sws.all.Rnomod-r_runpre.sws.all.Rnomod*r_prepost_3.sws.all.Rnomod)/sqrt((1-r_runpre.sws.all.Rnomod^2)*(1-r_prepost_3.sws.all.Rnomod^2)))^2;
	REV.post3.sws.all.Rnomod=((r_runpre.sws.all.Rnomod-r_runpost_3.sws.all.Rnomod*r_prepost_3.sws.all.Rnomod)/sqrt((1-r_runpost_3.sws.all.Rnomod^2)*(1-r_prepost_3.sws.all.Rnomod^2)))^2;
      else
	r_runpost_3.sws.all.Rmod=NaN;
	r_prepost_3.sws.all.Rmod=NaN;
	EV.post3.sws.all.Rmod=NaN;
	REV.post3.sws.all.Rmod=NaN;
	r_runpost_3.sws.all.Rnomod=NaN;
	r_prepost_3.sws.all.Rnomod=NaN;
	EV.post3.sws.all.Rnomod=NaN;
	REV.post3.sws.all.Rnomod=NaN;
      end
     
  %%% Pyr Only
    %%%%%% All ripple mod types
    r_runpost.sws.pyr.Rall=corrcoef(CorrMatrix.post.sws.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost.sws.pyr.Rall=r_runpost.sws.pyr.Rall(1,2);
    r_runpre.sws.pyr.Rall=corrcoef(CorrMatrix.pre.sws.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpre.sws.pyr.Rall=r_runpre.sws.pyr.Rall(1,2);
    r_prepost.sws.pyr.Rall=corrcoef(CorrMatrix.pre.sws.pyr.Rall,CorrMatrix.post.sws.pyr.Rall,'rows','complete');r_prepost.sws.pyr.Rall=r_prepost.sws.pyr.Rall(1,2);
    EV.sws.pyr.Rall=((r_runpost.sws.pyr.Rall-r_runpre.sws.pyr.Rall*r_prepost.sws.pyr.Rall)/sqrt((1-r_runpre.sws.pyr.Rall^2)*(1-r_prepost.sws.pyr.Rall^2)))^2;
    REV.sws.pyr.Rall=((r_runpre.sws.pyr.Rall-r_runpost.sws.pyr.Rall*r_prepost.sws.pyr.Rall)/sqrt((1-r_runpost.sws.pyr.Rall^2)*(1-r_prepost.sws.pyr.Rall^2)))^2;

    r_runpost_1.sws.pyr.Rall=corrcoef(CorrMatrix.post_1.sws.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost_1.sws.pyr.Rall=r_runpost_1.sws.pyr.Rall(1,2);
    r_runpost_2.sws.pyr.Rall=corrcoef(CorrMatrix.post_2.sws.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost_2.sws.pyr.Rall=r_runpost_2.sws.pyr.Rall(1,2);  
    r_prepost_1.sws.pyr.Rall=corrcoef(CorrMatrix.pre.sws.pyr.Rall,CorrMatrix.post_1.sws.pyr.Rall,'rows','complete');r_prepost_1.sws.pyr.Rall=r_prepost_1.sws.pyr.Rall(1,2);
    r_prepost_2.sws.pyr.Rall=corrcoef(CorrMatrix.pre.sws.pyr.Rall,CorrMatrix.post_2.sws.pyr.Rall,'rows','complete');r_prepost_2.sws.pyr.Rall=r_prepost_2.sws.pyr.Rall(1,2);   
    EV.post1.sws.pyr.Rall=((r_runpost_1.sws.pyr.Rall-r_runpre.sws.pyr.Rall*r_prepost_1.sws.pyr.Rall)/sqrt((1-r_runpre.sws.pyr.Rall^2)*(1-r_prepost_1.sws.pyr.Rall^2)))^2;
    REV.post1.sws.pyr.Rall=((r_runpre.sws.pyr.Rall-r_runpost_1.sws.pyr.Rall*r_prepost_1.sws.pyr.Rall)/sqrt((1-r_runpost_1.sws.pyr.Rall^2)*(1-r_prepost_1.sws.pyr.Rall^2)))^2;
    EV.post2.sws.pyr.Rall=((r_runpost_2.sws.pyr.Rall-r_runpre.sws.pyr.Rall*r_prepost_2.sws.pyr.Rall)/sqrt((1-r_runpre.sws.pyr.Rall^2)*(1-r_prepost_2.sws.pyr.Rall^2)))^2;
    REV.post2.sws.pyr.Rall=((r_runpre.sws.pyr.Rall-r_runpost_2.sws.pyr.Rall*r_prepost_2.sws.pyr.Rall)/sqrt((1-r_runpost_2.sws.pyr.Rall^2)*(1-r_prepost_2.sws.pyr.Rall^2)))^2;
    if ~isempty(CorrMatrix.post_3.sws.all.Rall)
      r_runpost_3.sws.pyr.Rall=corrcoef(CorrMatrix.post_3.sws.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost_3.sws.pyr.Rall=r_runpost_3.sws.pyr.Rall(1,2);
      r_prepost_3.sws.pyr.Rall=corrcoef(CorrMatrix.pre.sws.pyr.Rall,CorrMatrix.post_3.sws.pyr.Rall,'rows','complete');r_prepost_3.sws.pyr.Rall=r_prepost_3.sws.pyr.Rall(1,2);
      EV.post3.sws.pyr.Rall=((r_runpost_3.sws.pyr.Rall-r_runpre.sws.pyr.Rall*r_prepost_3.sws.pyr.Rall)/sqrt((1-r_runpre.sws.pyr.Rall^2)*(1-r_prepost_3.sws.pyr.Rall^2)))^2;
      REV.post3.sws.pyr.Rall=((r_runpre.sws.pyr.Rall-r_runpost_3.sws.pyr.Rall*r_prepost_3.sws.pyr.Rall)/sqrt((1-r_runpost_3.sws.pyr.Rall^2)*(1-r_prepost_3.sws.pyr.Rall^2)))^2;
    else
      r_runpost_3.sws.pyr.Rall=NaN;
      r_runpost_3.sws.pyr.Rall=NaN;
      EV.post3.sws.pyr.Rall=NaN;
      REV.post3.sws.pyr.Rall=NaN;
    end
    
    %%%% Ripple mod vs nonripple mod
    r_runpost.sws.pyr.Rmod=corrcoef(CorrMatrix.post.sws.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost.sws.pyr.Rmod=r_runpost.sws.pyr.Rmod(1,2);
    r_runpre.sws.pyr.Rmod=corrcoef(CorrMatrix.pre.sws.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpre.sws.pyr.Rmod=r_runpre.sws.pyr.Rmod(1,2);
    r_prepost.sws.pyr.Rmod=corrcoef(CorrMatrix.pre.sws.pyr.Rmod,CorrMatrix.post.sws.pyr.Rmod,'rows','complete');r_prepost.sws.pyr.Rmod=r_prepost.sws.pyr.Rmod(1,2);
    EV.sws.pyr.Rmod=((r_runpost.sws.pyr.Rmod-r_runpre.sws.pyr.Rmod*r_prepost.sws.pyr.Rmod)/sqrt((1-r_runpre.sws.pyr.Rmod^2)*(1-r_prepost.sws.pyr.Rmod^2)))^2;
    REV.sws.pyr.Rmod=((r_runpre.sws.pyr.Rmod-r_runpost.sws.pyr.Rmod*r_prepost.sws.pyr.Rmod)/sqrt((1-r_runpost.sws.pyr.Rmod^2)*(1-r_prepost.sws.pyr.Rmod^2)))^2;

    r_runpost_1.sws.pyr.Rmod=corrcoef(CorrMatrix.post_1.sws.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost_1.sws.pyr.Rmod=r_runpost_1.sws.pyr.Rmod(1,2);
    r_runpost_2.sws.pyr.Rmod=corrcoef(CorrMatrix.post_2.sws.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost_2.sws.pyr.Rmod=r_runpost_2.sws.pyr.Rmod(1,2);
    
    r_prepost_1.sws.pyr.Rmod=corrcoef(CorrMatrix.pre.sws.pyr.Rmod,CorrMatrix.post_1.sws.pyr.Rmod,'rows','complete');r_prepost_1.sws.pyr.Rmod=r_prepost_1.sws.pyr.Rmod(1,2);
    r_prepost_2.sws.pyr.Rmod=corrcoef(CorrMatrix.pre.sws.pyr.Rmod,CorrMatrix.post_2.sws.pyr.Rmod,'rows','complete');r_prepost_2.sws.pyr.Rmod=r_prepost_2.sws.pyr.Rmod(1,2);
    
    EV.post1.sws.pyr.Rmod=((r_runpost_1.sws.pyr.Rmod-r_runpre.sws.pyr.Rmod*r_prepost_1.sws.pyr.Rmod)/sqrt((1-r_runpre.sws.pyr.Rmod^2)*(1-r_prepost_1.sws.pyr.Rmod^2)))^2;
    REV.post1.sws.pyr.Rmod=((r_runpre.sws.pyr.Rmod-r_runpost_1.sws.pyr.Rmod*r_prepost_1.sws.pyr.Rmod)/sqrt((1-r_runpost_1.sws.pyr.Rmod^2)*(1-r_prepost_1.sws.pyr.Rmod^2)))^2;
    EV.post2.sws.pyr.Rmod=((r_runpost_2.sws.pyr.Rmod-r_runpre.sws.pyr.Rmod*r_prepost_2.sws.pyr.Rmod)/sqrt((1-r_runpre.sws.pyr.Rmod^2)*(1-r_prepost_2.sws.pyr.Rmod^2)))^2;
    REV.post2.sws.pyr.Rmod=((r_runpre.sws.pyr.Rmod-r_runpost_2.sws.pyr.Rmod*r_prepost_2.sws.pyr.Rmod)/sqrt((1-r_runpost_2.sws.pyr.Rmod^2)*(1-r_prepost_2.sws.pyr.Rmod^2)))^2;
    
    %%% Ripple nomod all
    r_runpost.sws.pyr.Rnomod=corrcoef(CorrMatrix.post.sws.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost.sws.pyr.Rnomod=r_runpost.sws.pyr.Rnomod(1,2);
    r_runpre.sws.pyr.Rnomod=corrcoef(CorrMatrix.pre.sws.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpre.sws.pyr.Rnomod=r_runpre.sws.pyr.Rnomod(1,2);
    r_prepost.sws.pyr.Rnomod=corrcoef(CorrMatrix.pre.sws.pyr.Rnomod,CorrMatrix.post.sws.pyr.Rnomod,'rows','complete');r_prepost.sws.pyr.Rnomod=r_prepost.sws.pyr.Rnomod(1,2);
    EV.sws.pyr.Rnomod=((r_runpost.sws.pyr.Rnomod-r_runpre.sws.pyr.Rnomod*r_prepost.sws.pyr.Rnomod)/sqrt((1-r_runpre.sws.pyr.Rnomod^2)*(1-r_prepost.sws.pyr.Rnomod^2)))^2;
    REV.sws.pyr.Rnomod=((r_runpre.sws.pyr.Rnomod-r_runpost.sws.pyr.Rnomod*r_prepost.sws.pyr.Rnomod)/sqrt((1-r_runpost.sws.pyr.Rnomod^2)*(1-r_prepost.sws.pyr.Rnomod^2)))^2;
    r_runpost_1.sws.pyr.Rnomod=corrcoef(CorrMatrix.post_1.sws.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost_1.sws.pyr.Rnomod=r_runpost_1.sws.pyr.Rnomod(1,2);
    r_runpost_2.sws.pyr.Rnomod=corrcoef(CorrMatrix.post_2.sws.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost_2.sws.pyr.Rnomod=r_runpost_2.sws.pyr.Rnomod(1,2);      
    r_prepost_1.sws.pyr.Rnomod=corrcoef(CorrMatrix.pre.sws.pyr.Rnomod,CorrMatrix.post_1.sws.pyr.Rnomod,'rows','complete');r_prepost_1.sws.pyr.Rnomod=r_prepost_1.sws.pyr.Rnomod(1,2);
    r_prepost_2.sws.pyr.Rnomod=corrcoef(CorrMatrix.pre.sws.pyr.Rnomod,CorrMatrix.post_2.sws.pyr.Rnomod,'rows','complete');r_prepost_2.sws.pyr.Rnomod=r_prepost_2.sws.pyr.Rnomod(1,2);    
    EV.post1.sws.pyr.Rnomod=((r_runpost_1.sws.pyr.Rnomod-r_runpre.sws.pyr.Rnomod*r_prepost_1.sws.pyr.Rnomod)/sqrt((1-r_runpre.sws.pyr.Rnomod^2)*(1-r_prepost_1.sws.pyr.Rnomod^2)))^2;
    REV.post1.sws.pyr.Rnomod=((r_runpre.sws.pyr.Rnomod-r_runpost_1.sws.pyr.Rnomod*r_prepost_1.sws.pyr.Rnomod)/sqrt((1-r_runpost_1.sws.pyr.Rnomod^2)*(1-r_prepost_1.sws.pyr.Rnomod^2)))^2;
    EV.post2.sws.pyr.Rnomod=((r_runpost_2.sws.pyr.Rnomod-r_runpre.sws.pyr.Rnomod*r_prepost_2.sws.pyr.Rnomod)/sqrt((1-r_runpre.sws.pyr.Rnomod^2)*(1-r_prepost_2.sws.pyr.Rnomod^2)))^2;
    REV.post2.sws.pyr.Rnomod=((r_runpre.sws.pyr.Rnomod-r_runpost_2.sws.pyr.Rnomod*r_prepost_2.sws.pyr.Rnomod)/sqrt((1-r_runpost_2.sws.pyr.Rnomod^2)*(1-r_prepost_2.sws.pyr.Rnomod^2)))^2;

    if ~isempty(CorrMatrix.post_3.sws.pyr.Rmod)
      r_runpost_3.sws.pyr.Rmod=corrcoef(CorrMatrix.post_3.sws.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost_3.sws.pyr.Rmod=r_runpost_3.sws.pyr.Rmod(1,2);
      r_prepost_3.sws.pyr.Rmod=corrcoef(CorrMatrix.pre.sws.pyr.Rmod,CorrMatrix.post_3.sws.pyr.Rmod,'rows','complete');r_prepost_3.sws.pyr.Rmod=r_prepost_3.sws.pyr.Rmod(1,2);
      EV.post3.sws.pyr.Rmod=((r_runpost_3.sws.pyr.Rmod-r_runpre.sws.pyr.Rmod*r_prepost_3.sws.pyr.Rmod)/sqrt((1-r_runpre.sws.pyr.Rmod^2)*(1-r_prepost_3.sws.pyr.Rmod^2)))^2;
      REV.post3.sws.pyr.Rmod=((r_runpre.sws.pyr.Rmod-r_runpost_3.sws.pyr.Rmod*r_prepost_3.sws.pyr.Rmod)/sqrt((1-r_runpost_3.sws.pyr.Rmod^2)*(1-r_prepost_3.sws.pyr.Rmod^2)))^2;
      r_runpost_3.sws.pyr.Rnomod=corrcoef(CorrMatrix.post_3.sws.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost_3.sws.pyr.Rnomod=r_runpost_3.sws.pyr.Rnomod(1,2);
      r_prepost_3.sws.pyr.Rnomod=corrcoef(CorrMatrix.pre.sws.pyr.Rnomod,CorrMatrix.post_3.sws.pyr.Rnomod,'rows','complete');r_prepost_3.sws.pyr.Rnomod=r_prepost_3.sws.pyr.Rnomod(1,2);
      EV.post3.sws.pyr.Rnomod=((r_runpost_3.sws.pyr.Rnomod-r_runpre.sws.pyr.Rnomod*r_prepost_3.sws.pyr.Rnomod)/sqrt((1-r_runpre.sws.pyr.Rnomod^2)*(1-r_prepost_3.sws.pyr.Rnomod^2)))^2;
      REV.post3.sws.pyr.Rnomod=((r_runpre.sws.pyr.Rnomod-r_runpost_3.sws.pyr.Rnomod*r_prepost_3.sws.pyr.Rnomod)/sqrt((1-r_runpost_3.sws.pyr.Rnomod^2)*(1-r_prepost_3.sws.pyr.Rnomod^2)))^2;
    else
      r_runpost_3.sws.pyr.Rmod=NaN;
      r_prepost_3.sws.pyr.Rmod=NaN;
      EV.post3.sws.pyr.Rmod=NaN;
      REV.post3.sws.pyr.Rmod=NaN;
      r_runpost_3.sws.pyr.Rnomod=NaN;
      r_prepost_3.sws.pyr.Rnomod=NaN;
      EV.post3.sws.pyr.Rnomod=NaN;
      REV.post3.sws.pyr.Rnomod=NaN;
    end
    
    
%%%%%%%%%%%% Ripple In/Out %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% All Cell Types  
    %%% All Rip Mod types
    r_runpost.rip_in.all.Rall=corrcoef(CorrMatrix.post.rip_in.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost.rip_in.all.Rall=r_runpost.rip_in.all.Rall(1,2);
    r_runpre.rip_in.all.Rall=corrcoef(CorrMatrix.pre.rip_in.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpre.rip_in.all.Rall=r_runpre.rip_in.all.Rall(1,2);
    r_prepost.rip_in.all.Rall=corrcoef(CorrMatrix.pre.rip_in.all.Rall,CorrMatrix.post.rip_in.all.Rall,'rows','complete');r_prepost.rip_in.all.Rall=r_prepost.rip_in.all.Rall(1,2);
    EV.rip_in.all.Rall=((r_runpost.rip_in.all.Rall-r_runpre.rip_in.all.Rall*r_prepost.rip_in.all.Rall)/sqrt((1-r_runpre.rip_in.all.Rall^2)*(1-r_prepost.rip_in.all.Rall^2)))^2;
    REV.rip_in.all.Rall=((r_runpre.rip_in.all.Rall-r_runpost.rip_in.all.Rall*r_prepost.rip_in.all.Rall)/sqrt((1-r_runpost.rip_in.all.Rall^2)*(1-r_prepost.rip_in.all.Rall^2)))^2;

    r_runpost.rip_out.all.Rall=corrcoef(CorrMatrix.post.rip_out.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpost.rip_out.all.Rall=r_runpost.rip_out.all.Rall(1,2);
    r_runpre.rip_out.all.Rall=corrcoef(CorrMatrix.pre.rip_out.all.Rall,CorrMatrix.run.all.Rall,'rows','complete');r_runpre.rip_out.all.Rall=r_runpre.rip_out.all.Rall(1,2);
    r_prepost.rip_out.all.Rall=corrcoef(CorrMatrix.pre.rip_out.all.Rall,CorrMatrix.post.rip_out.all.Rall,'rows','complete');r_prepost.rip_out.all.Rall=r_prepost.rip_out.all.Rall(1,2);
    EV.rip_out.all.Rall=((r_runpost.rip_out.all.Rall-r_runpre.rip_out.all.Rall*r_prepost.rip_out.all.Rall)/sqrt((1-r_runpre.rip_out.all.Rall^2)*(1-r_prepost.rip_out.all.Rall^2)))^2;
    REV.rip_out.all.Rall=((r_runpre.rip_out.all.Rall-r_runpost.rip_out.all.Rall*r_prepost.rip_out.all.Rall)/sqrt((1-r_runpost.rip_out.all.Rall^2)*(1-r_prepost.rip_out.all.Rall^2)))^2;
    
    %%%% Ripple mod vs nonripple mod
      % Mod
      r_runpost.rip_in.all.Rmod=corrcoef(CorrMatrix.post.rip_in.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost.rip_in.all.Rmod=r_runpost.rip_in.all.Rmod(1,2);
      r_runpre.rip_in.all.Rmod=corrcoef(CorrMatrix.pre.rip_in.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpre.rip_in.all.Rmod=r_runpre.rip_in.all.Rmod(1,2);
      r_prepost.rip_in.all.Rmod=corrcoef(CorrMatrix.pre.rip_in.all.Rmod,CorrMatrix.post.rip_in.all.Rmod,'rows','complete');r_prepost.rip_in.all.Rmod=r_prepost.rip_in.all.Rmod(1,2);
      EV.rip_in.all.Rmod=((r_runpost.rip_in.all.Rmod-r_runpre.rip_in.all.Rmod*r_prepost.rip_in.all.Rmod)/sqrt((1-r_runpre.rip_in.all.Rmod^2)*(1-r_prepost.rip_in.all.Rmod^2)))^2;
      REV.rip_in.all.Rmod=((r_runpre.rip_in.all.Rmod-r_runpost.rip_in.all.Rmod*r_prepost.rip_in.all.Rmod)/sqrt((1-r_runpost.rip_in.all.Rmod^2)*(1-r_prepost.rip_in.all.Rmod^2)))^2;

      r_runpost.rip_out.all.Rmod=corrcoef(CorrMatrix.post.rip_out.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpost.rip_out.all.Rmod=r_runpost.rip_out.all.Rmod(1,2);
      r_runpre.rip_out.all.Rmod=corrcoef(CorrMatrix.pre.rip_out.all.Rmod,CorrMatrix.run.all.Rmod,'rows','complete');r_runpre.rip_out.all.Rmod=r_runpre.rip_out.all.Rmod(1,2);
      r_prepost.rip_out.all.Rmod=corrcoef(CorrMatrix.pre.rip_out.all.Rmod,CorrMatrix.post.rip_out.all.Rmod,'rows','complete');r_prepost.rip_out.all.Rmod=r_prepost.rip_out.all.Rmod(1,2);
      EV.rip_out.all.Rmod=((r_runpost.rip_out.all.Rmod-r_runpre.rip_out.all.Rmod*r_prepost.rip_out.all.Rmod)/sqrt((1-r_runpre.rip_out.all.Rmod^2)*(1-r_prepost.rip_out.all.Rmod^2)))^2;
      REV.rip_out.all.Rmod=((r_runpre.rip_out.all.Rmod-r_runpost.rip_out.all.Rmod*r_prepost.rip_out.all.Rmod)/sqrt((1-r_runpost.rip_out.all.Rmod^2)*(1-r_prepost.rip_out.all.Rmod^2)))^2;   
      % No Mod
      r_runpost.rip_in.all.Rnomod=corrcoef(CorrMatrix.post.rip_in.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost.rip_in.all.Rnomod=r_runpost.rip_in.all.Rnomod(1,2);
      r_runpre.rip_in.all.Rnomod=corrcoef(CorrMatrix.pre.rip_in.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpre.rip_in.all.Rnomod=r_runpre.rip_in.all.Rnomod(1,2);
      r_prepost.rip_in.all.Rnomod=corrcoef(CorrMatrix.pre.rip_in.all.Rnomod,CorrMatrix.post.rip_in.all.Rnomod,'rows','complete');r_prepost.rip_in.all.Rnomod=r_prepost.rip_in.all.Rnomod(1,2);
      EV.rip_in.all.Rnomod=((r_runpost.rip_in.all.Rnomod-r_runpre.rip_in.all.Rnomod*r_prepost.rip_in.all.Rnomod)/sqrt((1-r_runpre.rip_in.all.Rnomod^2)*(1-r_prepost.rip_in.all.Rnomod^2)))^2;
      REV.rip_in.all.Rnomod=((r_runpre.rip_in.all.Rnomod-r_runpost.rip_in.all.Rnomod*r_prepost.rip_in.all.Rnomod)/sqrt((1-r_runpost.rip_in.all.Rnomod^2)*(1-r_prepost.rip_in.all.Rnomod^2)))^2;

      r_runpost.rip_out.all.Rnomod=corrcoef(CorrMatrix.post.rip_out.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpost.rip_out.all.Rnomod=r_runpost.rip_out.all.Rnomod(1,2);
      r_runpre.rip_out.all.Rnomod=corrcoef(CorrMatrix.pre.rip_out.all.Rnomod,CorrMatrix.run.all.Rnomod,'rows','complete');r_runpre.rip_out.all.Rnomod=r_runpre.rip_out.all.Rnomod(1,2);
      r_prepost.rip_out.all.Rnomod=corrcoef(CorrMatrix.pre.rip_out.all.Rnomod,CorrMatrix.post.rip_out.all.Rnomod,'rows','complete');r_prepost.rip_out.all.Rnomod=r_prepost.rip_out.all.Rnomod(1,2);
      EV.rip_out.all.Rnomod=((r_runpost.rip_out.all.Rnomod-r_runpre.rip_out.all.Rnomod*r_prepost.rip_out.all.Rnomod)/sqrt((1-r_runpre.rip_out.all.Rnomod^2)*(1-r_prepost.rip_out.all.Rnomod^2)))^2;
      REV.rip_out.all.Rnomod=((r_runpre.rip_out.all.Rnomod-r_runpost.rip_out.all.Rnomod*r_prepost.rip_out.all.Rnomod)/sqrt((1-r_runpost.rip_out.all.Rnomod^2)*(1-r_prepost.rip_out.all.Rnomod^2)))^2;   
      
  %%%%%%% Pyr only
    %%% All rip mod types
    r_runpost.rip_in.pyr.Rall=corrcoef(CorrMatrix.post.rip_in.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost.rip_in.pyr.Rall=r_runpost.rip_in.pyr.Rall(1,2);
    r_runpre.rip_in.pyr.Rall=corrcoef(CorrMatrix.pre.rip_in.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpre.rip_in.pyr.Rall=r_runpre.rip_in.pyr.Rall(1,2);
    r_prepost.rip_in.pyr.Rall=corrcoef(CorrMatrix.pre.rip_in.pyr.Rall,CorrMatrix.post.rip_in.pyr.Rall,'rows','complete');r_prepost.rip_in.pyr.Rall=r_prepost.rip_in.pyr.Rall(1,2);
    EV.rip_in.pyr.Rall=((r_runpost.rip_in.pyr.Rall-r_runpre.rip_in.pyr.Rall*r_prepost.rip_in.pyr.Rall)/sqrt((1-r_runpre.rip_in.pyr.Rall^2)*(1-r_prepost.rip_in.pyr.Rall^2)))^2;
    REV.rip_in.pyr.Rall=((r_runpre.rip_in.pyr.Rall-r_runpost.rip_in.pyr.Rall*r_prepost.rip_in.pyr.Rall)/sqrt((1-r_runpost.rip_in.pyr.Rall^2)*(1-r_prepost.rip_in.pyr.Rall^2)))^2;

    r_runpost.rip_out.pyr.Rall=corrcoef(CorrMatrix.post.rip_out.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpost.rip_out.pyr.Rall=r_runpost.rip_out.pyr.Rall(1,2);
    r_runpre.rip_out.pyr.Rall=corrcoef(CorrMatrix.pre.rip_out.pyr.Rall,CorrMatrix.run.pyr.Rall,'rows','complete');r_runpre.rip_out.pyr.Rall=r_runpre.rip_out.pyr.Rall(1,2);
    r_prepost.rip_out.pyr.Rall=corrcoef(CorrMatrix.pre.rip_out.pyr.Rall,CorrMatrix.post.rip_out.pyr.Rall,'rows','complete');r_prepost.rip_out.pyr.Rall=r_prepost.rip_out.pyr.Rall(1,2);
    EV.rip_out.pyr.Rall=((r_runpost.rip_out.pyr.Rall-r_runpre.rip_out.pyr.Rall*r_prepost.rip_out.pyr.Rall)/sqrt((1-r_runpre.rip_out.pyr.Rall^2)*(1-r_prepost.rip_out.pyr.Rall^2)))^2;
    REV.rip_out.pyr.Rall=((r_runpre.rip_out.pyr.Rall-r_runpost.rip_out.pyr.Rall*r_prepost.rip_out.pyr.Rall)/sqrt((1-r_runpost.rip_out.pyr.Rall^2)*(1-r_prepost.rip_out.pyr.Rall^2)))^2;

    %%%% Ripple mod vs nonripple mod
      % Mod
      r_runpost.rip_in.pyr.Rmod=corrcoef(CorrMatrix.post.rip_in.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost.rip_in.pyr.Rmod=r_runpost.rip_in.pyr.Rmod(1,2);
      r_runpre.rip_in.pyr.Rmod=corrcoef(CorrMatrix.pre.rip_in.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpre.rip_in.pyr.Rmod=r_runpre.rip_in.pyr.Rmod(1,2);
      r_prepost.rip_in.pyr.Rmod=corrcoef(CorrMatrix.pre.rip_in.pyr.Rmod,CorrMatrix.post.rip_in.pyr.Rmod,'rows','complete');r_prepost.rip_in.pyr.Rmod=r_prepost.rip_in.pyr.Rmod(1,2);
      EV.rip_in.pyr.Rmod=((r_runpost.rip_in.pyr.Rmod-r_runpre.rip_in.pyr.Rmod*r_prepost.rip_in.pyr.Rmod)/sqrt((1-r_runpre.rip_in.pyr.Rmod^2)*(1-r_prepost.rip_in.pyr.Rmod^2)))^2;
      REV.rip_in.pyr.Rmod=((r_runpre.rip_in.pyr.Rmod-r_runpost.rip_in.pyr.Rmod*r_prepost.rip_in.pyr.Rmod)/sqrt((1-r_runpost.rip_in.pyr.Rmod^2)*(1-r_prepost.rip_in.pyr.Rmod^2)))^2;

      r_runpost.rip_out.pyr.Rmod=corrcoef(CorrMatrix.post.rip_out.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpost.rip_out.pyr.Rmod=r_runpost.rip_out.pyr.Rmod(1,2);
      r_runpre.rip_out.pyr.Rmod=corrcoef(CorrMatrix.pre.rip_out.pyr.Rmod,CorrMatrix.run.pyr.Rmod,'rows','complete');r_runpre.rip_out.pyr.Rmod=r_runpre.rip_out.pyr.Rmod(1,2);
      r_prepost.rip_out.pyr.Rmod=corrcoef(CorrMatrix.pre.rip_out.pyr.Rmod,CorrMatrix.post.rip_out.pyr.Rmod,'rows','complete');r_prepost.rip_out.pyr.Rmod=r_prepost.rip_out.pyr.Rmod(1,2);
      EV.rip_out.pyr.Rmod=((r_runpost.rip_out.pyr.Rmod-r_runpre.rip_out.pyr.Rmod*r_prepost.rip_out.pyr.Rmod)/sqrt((1-r_runpre.rip_out.pyr.Rmod^2)*(1-r_prepost.rip_out.pyr.Rmod^2)))^2;
      REV.rip_out.pyr.Rmod=((r_runpre.rip_out.pyr.Rmod-r_runpost.rip_out.pyr.Rmod*r_prepost.rip_out.pyr.Rmod)/sqrt((1-r_runpost.rip_out.pyr.Rmod^2)*(1-r_prepost.rip_out.pyr.Rmod^2)))^2;   
      % No Mod
      r_runpost.rip_in.pyr.Rnomod=corrcoef(CorrMatrix.post.rip_in.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost.rip_in.pyr.Rnomod=r_runpost.rip_in.pyr.Rnomod(1,2);
      r_runpre.rip_in.pyr.Rnomod=corrcoef(CorrMatrix.pre.rip_in.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpre.rip_in.pyr.Rnomod=r_runpre.rip_in.pyr.Rnomod(1,2);
      r_prepost.rip_in.pyr.Rnomod=corrcoef(CorrMatrix.pre.rip_in.pyr.Rnomod,CorrMatrix.post.rip_in.pyr.Rnomod,'rows','complete');r_prepost.rip_in.pyr.Rnomod=r_prepost.rip_in.pyr.Rnomod(1,2);
      EV.rip_in.pyr.Rnomod=((r_runpost.rip_in.pyr.Rnomod-r_runpre.rip_in.pyr.Rnomod*r_prepost.rip_in.pyr.Rnomod)/sqrt((1-r_runpre.rip_in.pyr.Rnomod^2)*(1-r_prepost.rip_in.pyr.Rnomod^2)))^2;
      REV.rip_in.pyr.Rnomod=((r_runpre.rip_in.pyr.Rnomod-r_runpost.rip_in.pyr.Rnomod*r_prepost.rip_in.pyr.Rnomod)/sqrt((1-r_runpost.rip_in.pyr.Rnomod^2)*(1-r_prepost.rip_in.pyr.Rnomod^2)))^2;

      r_runpost.rip_out.pyr.Rnomod=corrcoef(CorrMatrix.post.rip_out.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpost.rip_out.pyr.Rnomod=r_runpost.rip_out.pyr.Rnomod(1,2);
      r_runpre.rip_out.pyr.Rnomod=corrcoef(CorrMatrix.pre.rip_out.pyr.Rnomod,CorrMatrix.run.pyr.Rnomod,'rows','complete');r_runpre.rip_out.pyr.Rnomod=r_runpre.rip_out.pyr.Rnomod(1,2);
      r_prepost.rip_out.pyr.Rnomod=corrcoef(CorrMatrix.pre.rip_out.pyr.Rnomod,CorrMatrix.post.rip_out.pyr.Rnomod,'rows','complete');r_prepost.rip_out.pyr.Rnomod=r_prepost.rip_out.pyr.Rnomod(1,2);
      EV.rip_out.pyr.Rnomod=((r_runpost.rip_out.pyr.Rnomod-r_runpre.rip_out.pyr.Rnomod*r_prepost.rip_out.pyr.Rnomod)/sqrt((1-r_runpre.rip_out.pyr.Rnomod^2)*(1-r_prepost.rip_out.pyr.Rnomod^2)))^2;
      REV.rip_out.pyr.Rnomod=((r_runpre.rip_out.pyr.Rnomod-r_runpost.rip_out.pyr.Rnomod*r_prepost.rip_out.pyr.Rnomod)/sqrt((1-r_runpost.rip_out.pyr.Rnomod^2)*(1-r_prepost.rip_out.pyr.Rnomod^2)))^2;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if strcmp(savevar,'on')
  save([xml '-ExplainedVariance-' structure '.mat'],'EV','REV','binsize','CorrMatrix','PvalMatrix','STindex','HPindex','structure','duration','is','STtype','HPtype');
end
