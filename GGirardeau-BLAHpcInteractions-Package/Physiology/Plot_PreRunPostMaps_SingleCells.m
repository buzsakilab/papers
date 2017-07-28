function Plot_PreRunPostMaps_SingleCells(celllist)

% Plots Firing Curves/maps and PETH for prerun, run, postrun for a cell or list of cells.
% Data explo utility

%PlotAirpuffRewardMod_SingleCells - Plots Maps, Curves, PETH for a list of cells [rat session shank cellnumber]
%
%  USAGE
%
%    PlotAirpuffRewardMod_SingleCells (celllist,<options>)
%
%    celllist       rat/sess/sahnk/cell                              
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
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

load('/media/Data-01/All-Rats/sessionindexing.mat');


struc=input('Structure?','s')

sessionlist=unique(celllist(:,1:2),'rows');
for i=1:size(sessionlist,1)
  ratsess=sessionlist(i,:);
  path=xmlpath{ismember(ratsessionindex,ratsess,'rows')}
  cd(path)
  xml=path(end-14:end-1)
  sessioncells=celllist(celllist(:,1)==ratsess(1)&celllist(:,2)==ratsess(2),3:4)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mapping Parameters (output from AirpuffRewardMod)
  if exist([xml '-AirpuffRewardMod-2.mat'])==2
    load([xml '-AirpuffRewardMod-2.mat'])
    load([xml '-AirpuffRewardMod-SafeRew.mat'])
  else
    warning(['No Data for the session : ' xml]);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ids=find(ismember(Index(:,1:2),sessioncells,'rows'));
  
  for i=1:length(ids)
    id=ids(i);
    if ~isempty(Spiketimes{ids(i)}.run);
      [spikepos,~]=Interpolate(Pos.runpos,Spiketimes{ids(i)}.run);
    else
      spikepos=[];
    end
    figure('Position',[180 291 1690 454]);
    subplot(3,5,1:2)
    PlotColorMap(MazeMap{id}.prerun.map.rate,'newfig','off','ydir','reverse');
    ylim([0 60])
    colorbar 
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    subplot(3,5,6:7)
    PlotColorMap(MazeMap{id}.run.map.rate,'newfig','off','ydir','reverse');
    ylim([0 60])
    colorbar
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
    subplot(3,5,11:12)
    PlotColorMap(MazeMap{id}.postrun.map.rate,'newfig','off','ydir','reverse');
    ylim([0 60])
    colorbar 
    PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));

    subplot(3,5,3:4)
    bar(MazeCurve{id}.prerun.curve.x,MazeCurve{id}.prerun.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
    xlim([0 1])
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    subplot(3,5,8:9)
    bar(MazeCurve{id}.run.curve.x,MazeCurve{id}.run.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
    xlim([0 1])
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    subplot(3,5,13:14)
    bar(MazeCurve{id}.postrun.curve.x,MazeCurve{id}.postrun.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
    xlim([0 1])
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    
    [m3r,t3r]=SyncHist(RewardPETH{id}.prerun.raster,RewardPETH{id}.prerun.indices,'mode','sum','smooth',0,'durations',[-10 10]);
  %    [m3r,t3r]=SyncHist(RewardPETH{id}.prerun.raster,RewardPETH{id}.prerun.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    subplot(3,5,5)
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
    end
    
  %    [m3r,t3r]=SyncHist(RewardPETH{id}.run.raster,RewardPETH{id}.run.indices,'mode','sum','smooth',0,'durations',[-10 10]);
    [m3r,t3r]=SyncHist(SafeRewPETH{id}.run.raster,SafeRewPETH{id}.run.indices,'mode','sum','smooth',0,'durations',[-10 10]);
  %    [m3r,t3r]=SyncHist(RewardPETH{id}.saferun.raster,RewardPETH{id}.saferun.indices,'mode','sum','smooth',0,'durations',[-10 10]);
    subplot(3,5,10)
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
    end 
    
    [m3r,t3r]=SyncHist(SafeRewPETH{id}.postrun.raster,SafeRewPETH{id}.postrun.indices,'mode','sum','smooth',0,'durations',[-10 10]);
  %    [m3r,t3r]=SyncHist(RewardPETH{id}.postrun.raster,RewardPETH{id}.postrun.indices,'mode','sum','smooth',0,'durations',[-10 10]);
  %    [m3r,t3r]=SyncHist(RewardPETH{id}.postrun.raster,RewardPETH{id}.postrun.indices,'mode','sum','smooth',0,'durations',[-2 2]);
    subplot(3,5,15)
    
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
    end
    
    suptitle([xml ' - cell ' int2str(sessioncells(i,1)) ' shank ' int2str(sessioncells(i,2))]);
    
    
    cd('/media/Data-01/All-Rats/AllRats-EVcontributingCells')    
    saveas(gcf,[struc ' - ' xml ' - Shank ' int2str(sessioncells(i,2)) ' Cell ' int2str(sessioncells(i,1))],'png');
    close
  end
end










