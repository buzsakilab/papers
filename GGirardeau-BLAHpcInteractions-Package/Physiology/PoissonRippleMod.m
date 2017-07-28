function poissonripplemod = PoissonRippleMod(session,varargin)

%PoissonRippleMod - Calculate the significance of global firing rate modulation of individual cells by ripples using a Poisson Test. Plots PETH for ripples (1 fig/shank)
%
%  USAGE
%
%    poissonripplemod = PoissonRippleMod (session,<options>)
%
%    session		recording session
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'DB'		DataBase name if figures/variables are to be stored in DB (default = [])
%     'savevar'		Save matlab variable (Default 'off')
%     'bigripples'      Use only large ripples (Default 'off')
%     'savepng'         Save png figure default 'off'
%    =========================================================================
%
%  OUTPUT
%
%    poissonripplemod : matrix [Shank / Cell / ID / pIncrease / pDecrease / Surprise]
%
%  NOTES
%
%  SEE
%
%    See also : Eran's PoissonTest
%
% Jan 2014 by Gabrielle Girardeau calling Eran Stark's PoissonTest.m
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
DB=[];
savevar ='off';
bigripples='off';
savepng = 'off';
% For plotting only
pval=0.001;
duration = 2;
binSize = 0.01;
smooth=2;


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
		case 'db'
			DB = varargin{i+1};
		case 'show'
		       show = varargin{i+1};
		case 'bigripples'
		       bigripples = varargin{i+1};
		case 'savevar'
		       savevar = varargin{i+1};
		case 'pval'
		       pval = varargin{i+1};
		case 'duration'
		      duration = varargin{i+1};
		case 'smooth'
		      smooth = varargin{i+1};
		case 'binSize'
		      binSize = varargin{i+1};
		case 'savepng'
		      savepng = varargin{i+1};
		otherwise,
			error(['Unknown property ''' num2str(varargin{i})]);
	end
end

cd(session);
xml=session(end-13:end);
SetCurrentSession(xml);

ripples=GetRippleEvents;
channel=GetRippleChannel;

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
index=unique(spikes(:,2:4),'rows'); % shank / cell/ ID

load([xml '-States.mat']);

% Restrict to sws
ripples=Restrict(ripples,sws);
spikes=Restrict(spikes,sws);

% baseline epochs
bufferedripples=[ripples(:,1)-0.1 ripples(:,3)+0.1];
[baseline,ind]=SubtractIntervals(sws,bufferedripples,'strict','on');
totalbaselinetime=sum(baseline(:,2)-baseline(:,1));

baselinespikes=Restrict(spikes,baseline);
nbaselinespikes=size(baselinespikes,1);

% Restrict further to the third of ripples with the highest amplitude
if strcmp(bigripples,'on')
  [maps,data]=RippleStats(fil,ripples);
  [sorted,indx]=sortrows(data.peakAmplitude);
  third=round(size(data.peakAmplitude,1)/3);
  idxhigh=indx(end-third:end);
  ripples=ripples(idxhigh,:);
  ripples=sort(ripples);
end

totalrippletime=sum(ripples(:,3)-ripples(:,1));
ripplespikes=Restrict(spikes,[ripples(:,1) ripples(:,3)]);
nripplespikes=size(ripplespikes,1);

for i=1:size(index,1)
  ncellbaselinespikes=length(baselinespikes(baselinespikes(:,4)==i,1));
  ncellripplespikes=length(ripplespikes(ripplespikes(:,4)==i,1));
  if ncellbaselinespikes~=0 & ncellripplespikes~=0
    [pInc(i) pDec(i) surp(i)] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
  else
    pInc(i)=NaN;
    pDec(i)=NaN;
    surp(i)=NaN;
  end
end

poissonripplemod=[index pInc' pDec' surp'];
if strcmp(savevar,'on')
  save([xml '-PoissonRippleMod.mat'],'poissonripplemod')
end

% Plot
figuresize=[453 85 1331 839];
rippleccg=[ripples(:,2) ones(size(ripples,1),1) ones(size(ripples,1),1)];
shanklist=unique(index(:,1));
for j = 1:length(shanklist)
  shank=shanklist(j);
  block = spikes(spikes(:,2)==shank,:);
  if ~isempty(block)
%        u=figure;
      figure;
%        set(u,'Position',figuresize);
      set(gcf,'Position',figuresize);
      suptitle([xml 'Shank ' int2str(shank)]);
      for cell=2:length(unique(block(:,3)))+1;
	  cellspiketimes=block(block(:,3)==cell,1);
	  if ~isempty(cellspiketimes)
	    spikesccg=[cellspiketimes 2*ones(length(cellspiketimes),1) 2*ones(length(cellspiketimes),1)];
	    forccg=[rippleccg; spikesccg];
	    forccg=sortrows(forccg,1);
	    [ccg,ti]=CCG(forccg(:,1),forccg(:,2),'binSize',binSize,'duration',duration,'smooth',smooth);
	    SquareSubplot(length(unique(block(:,3))),cell-1);
	    bar(ti,ccg(:,1,2));
	    PlotHVLines(0,'v',':r');
	    xlim([min(ti) max(ti)]);
	    xlabel(['Cell # ' int2str(cell)]);
	    if poissonripplemod(poissonripplemod(:,1)==shank&poissonripplemod(:,2)==cell,4)<pval/2
	      PlotColorIntervals([-1 1],'DarkSeaGreen','v');
	    elseif poissonripplemod(poissonripplemod(:,1)==shank&poissonripplemod(:,2)==cell,5)<pval/2
	      PlotColorIntervals([-1 1],'LightSalmon','v');
	    end
	  end
      end

      if (strcmp(savepng,'on'))
	print(gcf,['/media/Data-01/All-Rats/AllRats-PoissonRippleMod/images/' xml '-Shank' int2str(shank)],'-dpng');
      end

      % DataBase
      if ~ isempty(DB)
	DBConnect
	DBUse(DB);
	DBAddFigure(gcf,['Hpc-Amygdala'],[xml ' shank ' int2str(shank) '  - Poisson Ripple Modulation'],'','','PoissonRippleMod.m');
	close(gcf);
      end
   end
end








