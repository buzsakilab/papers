function PlotAirpuffRewardMod_SingleCells(celllist,varargin)

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

sessionlist=unique(celllist(:,1:2),'rows');

struc=input('Structure?','s')

for i=1:size(sessionlist,1)
  ratsess=sessionlist(i,:);
  path=xmlpath{ismember(ratsessionindex,ratsess,'rows')}
  cd(path)
  xml=path(end-14:end-1)
  sessioncells=celllist(celllist(:,1)==ratsess(1)&celllist(:,2)==ratsess(2),3:4)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mapping Parameters (output from AirpuffRewardMod)
  if exist([xml '-AirpuffRewardMod-2.mat'])==2
    load([xml '-AirpuffRewardMod-2.mat'])
  else
    warning(['No Data for the session : ' xml]);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ids=find(ismember(Index(:,1:2),sessioncells,'rows'));
  
  for i=1:length(ids)
    if ~isempty(Spiketimes{ids(i)}.run);
      [spikepos,~]=Interpolate(Pos.runpos,Spiketimes{ids(i)}.run);
    else
      spikepos=[];
    end

    fg=figure;
    set(fg,'Position',[15 60 1903 910]);
    subplot(5,7,[1 8 15 22 29])
    plot(Pos.runpos(:,2),Pos.runpos(:,1),'Color',rgb('Gray'));
    hold on
    if ~isempty(spikepos)
      plot(spikepos(:,2),spikepos(:,1),'.','Color',rgb('Crimson'));
    end
    PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
    
    %%% 2D maps
    subplot(5,7,2:3)
      PlotColorMap(MazeMap{ids(i)}.run.map.rate,'newfig','off','ydir','reverse');
      ylim([0 60]);
      PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
      colorbar;
      col=clim;
    subplot(5,7,9:10)
      PlotColorMap(MazeMap{ids(i)}.run.safelapsmap.rate,'newfig','off','ydir','reverse');
      ylim([0 60]);
      PlotHVLines(Airpuff.loc*137,'v','Color',rgb('Gray'));
      colorbar;
      set(gca,'clim',col);
    subplot(5,7,16:17)
      PlotColorMap(MazeMap{ids(i)}.run.airpufflapsmap.rate,'newfig','off','ydir','reverse');
      ylim([0 60]);
      PlotHVLines(Airpuff.loc*137,'v','Color',rgb('OrangeRed'));
      colorbar;
      set(gca,'clim',col);

    %%% 1D curves - Rates
    subplot(5,7,4:5)
      bar(MazeCurve{ids(i)}.run.curve.x,MazeCurve{ids(i)}.run.curve.rate,'FaceColor',rgb('SeaGreen'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
      xlabel(['RUN - Field peak: ' num2str(MazeCurve{ids(i)}.run.stats.peak) 'Hz. Field mean rate: ' num2str(MazeCurve{ids(i)}.run.stats.mean) 'Hz']);
      ylabel('Spike Rate')
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      if strcmp(Airpuff.dir,'LtoR')
	arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      elseif strcmp(Airpuff.dir,'RtoL')
	arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      end
    subplot(5,7,11:12)
      bar(MazeCurve{ids(i)}.run.safelapscurve.x,MazeCurve{ids(i)}.run.safelapscurve.rate,'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('Gray'));
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      ylabel('Spike Rate')
    subplot(5,7,18:19)
      bar(MazeCurve{ids(i)}.run.airpufflapscurve.x,MazeCurve{ids(i)}.run.airpufflapscurve.rate,'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
      ylabel('Spike Rate')
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      if strcmp(Airpuff.dir,'LtoR')
	arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      elseif strcmp(Airpuff.dir,'RtoL')
	arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      end

      %%% 1D curves - SPike Counts
      subplot(5,7,6:7)
      bar(MazeCurve{ids(i)}.run.curve.x,MazeCurve{ids(i)}.run.curve.count,'FaceColor',rgb('DarkCyan'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
      xlabel(['RUN - Field peak: ' num2str(MazeCurve{ids(i)}.run.stats.peak) 'Hz. Field mean rate: ' num2str(MazeCurve{ids(i)}.run.stats.mean) 'Hz']);
      ylabel('Spike Count')
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      if strcmp(Airpuff.dir,'LtoR')
	arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      elseif strcmp(Airpuff.dir,'RtoL')
	arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      end
    subplot(5,7,13:14)
      bar(MazeCurve{ids(i)}.run.safelapscurve.x,MazeCurve{ids(i)}.run.safelapscurve.count,'FaceColor',rgb('DarkTurquoise'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('Gray'));
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      ylabel('Spike Count')
    subplot(5,7,20:21)
      bar(MazeCurve{ids(i)}.run.airpufflapscurve.x,MazeCurve{ids(i)}.run.airpufflapscurve.count,'FaceColor',rgb('DarkTurquoise'),'EdgeColor','none','BarWidth',1);
      hold on;
      xlim([0 1]);
      PlotHVLines(Airpuff.loc,'v','Color',rgb('OrangeRed'));
      ylabel('Spike Count')
      % Plot Airpuff direction arrow
      ylimit=get(gca,'YLim');
      if strcmp(Airpuff.dir,'LtoR')
	arrow([Airpuff.loc-0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      elseif strcmp(Airpuff.dir,'RtoL')
	arrow([Airpuff.loc+0.1 ylimit(2)*0.8],[Airpuff.loc ylimit(2)*0.8],'Length',5);
      end


      % Airpuff PETHs
      [m1a,t1a]=SyncHist(AirpuffPETH{ids(i)}.raster,AirpuffPETH{ids(i)}.indices,'mode','sum','smooth',0);
      [m2a,t2a]=SyncHist(AirpuffPETH{ids(i)}.raster,AirpuffPETH{ids(i)}.indices,'mode','sum','smooth',0,'durations',[-2 2]);
      [m3a,t3a]=SyncHist(AirpuffPETH{ids(i)}.raster,AirpuffPETH{ids(i)}.indices,'mode','sum','smooth',0,'durations',[-10 10]);

      % Reward PETHs
      % All runs
      [m1r,t1r]=SyncHist(RewardPETH{ids(i)}.allrun.raster,RewardPETH{ids(i)}.allrun.indices,'mode','sum','smooth',0);
      [m2r,t2r]=SyncHist(RewardPETH{ids(i)}.allrun.raster,RewardPETH{ids(i)}.allrun.indices,'mode','sum','smooth',0,'durations',[-2 2]);
      [m3r,t3r]=SyncHist(RewardPETH{ids(i)}.allrun.raster,RewardPETH{ids(i)}.allrun.indices,'mode','sum','smooth',0,'durations',[-10 10]);

    subplot(5,7,25)
    if ~isempty(m1a)
      bar(t1a,m1a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      ylabel('Airpuff PETH');
      xlabel('500 ms');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,7,26)
    if ~isempty(m2a)
      bar(t2a,m2a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','AirpuffMod');
    end
    subplot(5,7,27)
    if ~isempty(m3a)
      bar(t3a,m3a,'FaceColor',rgb('Crimson'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','AirpuffMod');
    end
    hold on;
    subplot(5,7,32)
    if ~isempty(m1a)
      bar(t1r,m1r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      ylabel('Reward PETH');
      xlabel('500 ms');
      set(gca,'Tag','RewardMod');
    end
    subplot(5,7,33)
    if ~isempty(m2a)
      bar(t2r,m2r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('2 sec');
      set(gca,'Tag','RewardMod');
    end
    subplot(5,7,34)
    if ~isempty(m3r)
      bar(t3r,m3r,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none','BarWidth',1);
      xlabel('10 sec');
      set(gca,'Tag','RewardMod');
    end
    

   
    suptitle([xml ' - Shank ' int2str(Index(Index(:,3)==ids(i),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i),2)) ' (ID: ' int2str(ids(i))]);
%      cd('/media/Data-01/All-Rats/AllRats-EVcontributingCells')
    
%      saveas(gcf,[struc ' - ' xml ' - Shank ' int2str(Index(Index(:,3)==ids(i),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i),2))],'png');
%      plot2svg([struc ' - ' xml ' - Shank ' int2str(Index(Index(:,3)==ids(i),1)) ' Cell ' int2str(Index(Index(:,3)==ids(i),2)) '.svg'],gcf);
%      close
  end
end












