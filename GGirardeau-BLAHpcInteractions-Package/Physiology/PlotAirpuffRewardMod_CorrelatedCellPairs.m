function PlotAirpuffRewardMod_CorrelatedCellPairs(xml,diffthreshold,runthreshold,varargin)

%PlotAirpuffRewardMod_CorrelatedCellPairs - Plots Maps, Curves, PETH for correlated hpc-structure cell pairs
%
%  USAGE
%
%    PlotAirpuffRewardMod_CorrelatedCellPairs (session,diffthreshold,runthreshold,<options>)
%
%    session        session, xml or full path                              
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%       diffthreshold	threshold for pre/post correlation difference to select cells to plot
%	runthreshold	threshold for run correlation value to select cells to plot
%	celllist	list of cell pairs to plot (doublets of indices) or qudruplet of shank/cell pair. (ignoring correlation selection criteria for slecteing cells to plot)
%	savefig		'on','off'(Default)
%    =========================================================================
%
%  OUTPUT
%
%    Figures
%
%  NOTE
%
%  SEE
%
%    See also : FiringCurve, MapStats (FMAToolbox)
%
% Gabrielle Girardeau, November 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
cellindices = [];
diffthreshold = 0.015;
runthreshold = [];

xml=[];

% Check number of inputs
if nargin < 2,
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
			cellindices = varargin{i+1};
		case 'diffthreshold',
			diffthreshold = varargin{i+1};
		case 'runthreshold',
			runthreshold = varargin{i+1};	
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mapping Parameters (output from AirpuffRewardMod)
if exist([xml '-AirpuffRewardMod-2.mat'])==2
  load([xml '-AirpuffRewardMod-2.mat'])
else
  warning(['No Data for the session : ' xml]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([xml '-ExplainedVariance-BLA.mat'])

DiffMatrix=CorrMatrix.post.sws.all.Rall-CorrMatrix.pre.sws.all.Rall;

if isempty(runthreshold) & ~isempty(diffthreshold)
  selection=DiffMatrix>diffthreshold;
elseif ~isempty(runthreshold) & isempty(diffthreshold)
  selection=CorrMatrix.run.all.Rall>runthreshold
elseif ~isempty(runthreshold) & ~isempty(diffthreshold)
  selection=DiffMatrix>diffthreshold & CorrMatrix.run.all.Rall>runthreshold;
end

if isempty(cellindices)
  [i,j]=find(selection);
  preCorr=CorrMatrix.pre.sws.all.Rall(selection);
  postCorr=CorrMatrix.post.sws.all.Rall(selection);
  runCorr=CorrMatrix.run.all.Rall(selection);
else
  i=cellindices(:,1)
  j=cellindices(:,2)
  for k=1:size(cellindices,1)
    preCorr(k)=CorrMatrix.pre.sws.all.Rall(i(k),j(k));
    postCorr(k)=CorrMatrix.post.sws.all.Rall(i(k),j(k));
    runCorr(k)=CorrMatrix.run.all.Rall(i(k),j(k));
  end
end
celllist=[HPindex(i,:) STindex(j,:)];
typelist=[HPtype(i,:) STtype(j,:)];
  
for i=1:size(celllist,1)
  ids(i,1)=find(ismember(Index(:,1:2),celllist(i,1:2),'rows'));
end
for i=1:size(celllist,1)
  ids(i,2)=find(ismember(Index(:,1:2),celllist(i,4:5),'rows'));
end
types=typelist;

for i=1:size(ids,1)
  if size((Spiketimes{ids(i,1)}.run),1)>0
    [spikepos1,~]=Interpolate(Pos.runpos,Spiketimes{ids(i,1)}.run);
    [spikepos2,~]=Interpolate(Pos.runpos,Spiketimes{ids(i,2)}.run);

    fg=figure;
    set(fg,'Position',[15 60 1903 910]);
    
    subplot(5,5,[1 6 11 16])
  %    plot(Pos.runpos(:,2),Pos.runpos(:,1),'Color',rgb('Gray'));
  %    hold on
  %    plot(spikepos1(:,2),spikepos1(:,1),'.','Color',rgb('Crimson'));
  %    plot(spikepos2(:,2),spikepos2(:,1),'.','Color',rgb('MidnightBlue'));
  %    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    
    subplot(5,5,21)
    plot([1 2 3],[preCorr(i) runCorr(i) postCorr(i)],'.-k','MarkerSize',25)
    xlim([0 4]);
    xlabel(['PrePost diff = ' num2str(postCorr(i)-preCorr(i))]);
    set(gca,'XtickLabel',[{'' 'Pre' 'Run' 'Post' ''}])
    
    %%% 2D maps
    subplot(5,5,2:3)
    PlotColorMap(MazeMap{ids(i,1)}.run.map.rate,'newfig','off','ydir','reverse');
    ylim([0 60]);
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    colorbar;
    col=clim;

    if types(i,1)==1
      title([xml ' - Shank ' int2str(Index(Index(:,3)==ids(i,1),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i,1),2)) ' (ID: ' int2str(ids(i,1)) ') - put. PYR']);
    elseif types(i,1)==2
      title([xml ' - Shank ' int2str(Index(Index(:,3)==ids(i,1),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i,1),2)) ' (ID: ' int2str(ids(i,1)) ') - put. INT']);
    end
    
    subplot(5,5,12:13)
    PlotColorMap(MazeMap{ids(i,1)}.run.airpufflapsmap.rate,'newfig','off','ydir','reverse');
    ylim([0 60]);
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    colorbar;
    set(gca,'clim',col);
    
    %%% 1D curves - Rates
    subplot(5,5,7:8)
    bar(MazeCurve{ids(i,1)}.run.curve.x,MazeCurve{ids(i,1)}.run.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
    hold on;
    xlim([0 1]);
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    xlabel(['RUN - Field peak: ' num2str(MazeCurve{ids(i,1)}.run.stats.peak) 'Hz. Field mean rate: ' num2str(MazeCurve{ids(i,1)}.run.stats.mean) 'Hz']);
    ylabel('Spike Rate')
    % Plot Airpuff direction arrow
    ylimit=get(gca,'YLim');
    if strcmp(Airpuff.dir,'LtoR')
      arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
    elseif strcmp(Airpuff.dir,'RtoL')
      arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
    end
    colorbar
    
    % Airpuff PETHs
    [m2a,t2a]=SyncHist(AirpuffPETH{ids(i,1)}.raster,AirpuffPETH{ids(i,1)}.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    [m3a,t3a]=SyncHist(AirpuffPETH{ids(i,1)}.raster,AirpuffPETH{ids(i,1)}.indices,'mode','sum','smooth',0,'durations',[-10 10]);
    % Reward PETHs
    % All runs
    [m2r,t2r]=SyncHist(RewardPETH{ids(i,1)}.allrun.raster,RewardPETH{ids(i,1)}.allrun.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    [m3r,t3r]=SyncHist(RewardPETH{ids(i,1)}.allrun.raster,RewardPETH{ids(i,1)}.allrun.indices,'mode','sum','smooth',0,'durations',[-10 10]);

    subplot(5,5,17)
    if ~isempty(m2a)
      bar(t2a,m2a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,5,18)
    if ~isempty(m3a)
      bar(t3a,m3a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,5,22)
    if ~isempty(m2a)
      bar(t2r,m2r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','RewardMod');
    end
    subplot(5,5,23)
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','RewardMod');
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second cell  
    %%% 2D maps
    subplot(5,5,4:5)
    PlotColorMap(MazeMap{ids(i,2)}.run.map.rate,'newfig','off','ydir','reverse');
    ylim([0 60]);
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    colorbar;
    col=clim;

    if types(i,2)==1
      title([xml ' - Shank ' int2str(Index(Index(:,3)==ids(i,2),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i,2),2)) ' (ID: ' int2str(ids(i,2)) ') - put. PYR']);
    elseif types(i,2)==2
      title([xml ' - Shank ' int2str(Index(Index(:,3)==ids(i,2),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i,2),2)) ' (ID: ' int2str(ids(i,2)) ') - put. INT']);
    end
    
    subplot(5,5,14:15)
    PlotColorMap(MazeMap{ids(i,2)}.run.airpufflapsmap.rate,'newfig','off','ydir','reverse');
    ylim([0 60]);
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    colorbar;
    set(gca,'clim',col);

    %%% 1D curves - Rates
    subplot(5,5,9:10)
    bar(MazeCurve{ids(i,2)}.run.curve.x,MazeCurve{ids(i,2)}.run.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
    hold on;
    xlim([0 1]);
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    xlabel(['RUN - Field peak: ' num2str(MazeCurve{ids(i,2)}.run.stats.peak) 'Hz. Field mean rate: ' num2str(MazeCurve{ids(i,2)}.run.stats.mean) 'Hz']);
    ylabel('Spike Rate')
    % Plot Airpuff direction arrow
    ylimit=get(gca,'YLim');
    if strcmp(Airpuff.dir,'LtoR')
      arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
    elseif strcmp(Airpuff.dir,'RtoL')
      arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
    end
    colorbar
    % Airpuff PETHs
    [m2a,t2a]=SyncHist(AirpuffPETH{ids(i,2)}.raster,AirpuffPETH{ids(i,2)}.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    [m3a,t3a]=SyncHist(AirpuffPETH{ids(i,2)}.raster,AirpuffPETH{ids(i,2)}.indices,'mode','sum','smooth',0,'durations',[-10 10]);
    % Reward PETHs
    % All runs
    [m2r,t2r]=SyncHist(RewardPETH{ids(i,2)}.allrun.raster,RewardPETH{ids(i,2)}.allrun.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    [m3r,t3r]=SyncHist(RewardPETH{ids(i,2)}.allrun.raster,RewardPETH{ids(i,2)}.allrun.indices,'mode','sum','smooth',0,'durations',[-10 10]);

    subplot(5,5,19)
    if ~isempty(m2a)
      bar(t2a,m2a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,5,20)
    if ~isempty(m3a)
      bar(t3a,m3a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,5,24)
    if ~isempty(m2a)
      bar(t2r,m2r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','RewardMod');
    end
    subplot(5,5,25)
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','RewardMod');
    end  
  end
end












