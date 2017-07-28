function [safemap,apmap,airP] = VelocityMaps(session)

%VelocityMaps - Calculates speed curves in the safe and airpuff directions for single session
%
%  USAGE
%
%    [safemap,apmap,airP] = VelocityMaps(session)
%
%    session    path to session
%
%  OUTPUT
%  
%   safemap         speed curves for safe runs (.prerun .run .postrun)
%   apmap           speed curves for airpuff runs
%   airP            normalized airpuff position
%
%  SEE
%
%    VelocityMaps_All
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


cd(session);
xml=session(end-13:end);
SetCurrentSession([session '/' xml],'spikes','off');

load([xml '-LapType.mat']);
%  prev=load('airpuffPrevious.mat');
load('runintervals.mat');
load('airpuff.mat');
pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);

vel=LinearVelocity(pos);
smvel=[vel(:,1) Smooth(vel(:,2),2)];

velpos=[pos(:,1:2) smvel(:,2)];

%  figure;
%  subplot(3,1,2);hold on;
%  for i=1:length(aplaps.run)
%      lapvelpos=Restrict(velpos,aplaps.run(i,:));
%      PlotXY(lapvelpos(:,2:3));
%  end
%  subplot(3,1,1);hold on;
%  for i=1:length(aplaps.prerun)
%      lapvelpos=Restrict(velpos,aplaps.prerun(i,:));
%      PlotXY(lapvelpos(:,2:3));
%  end
%  subplot(3,1,3);hold on;
%  for i=1:length(aplaps.postrun)
%      lapvelpos=Restrict(velpos,aplaps.postrun(i,:));
%      PlotXY(lapvelpos(:,2:3));
%  end


velposap.prerun=Restrict(velpos,aplaps.prerun);
velposap.run=Restrict(velpos,aplaps.run);
velpossafe.prerun=Restrict(velpos,safelaps.prerun);
velpossafe.run=Restrict(velpos,safelaps.run);

velposap.postrun.all=Restrict(velpos,aplaps.postrun);
velpossafe.postrun.all=Restrict(velpos,safelaps.postrun);

keyboard

if size(aplaps.postrun,1)>=5 && size(safelaps.postrun,1)>=5
    velposap.postrun.firsts=Restrict(velpos,aplaps.postrun(1:5,:));
    velpossafe.postrun.firsts=Restrict(velpos,safelaps.postrun(1:5,:));
else
    velposap.postrun.firsts=Restrict(velpos,aplaps.postrun);
    velpossafe.postrun.firsts=Restrict(velpos,safelaps.postrun);
end
    
%  [x,airP,prevairP] = ZeroToOne(velposap.run(:,2),airpuff.loc*0.43,prev.airpuff.loc*0.43);
[x,airP] = ZeroToOne(velposap.run(:,2),airpuff.loc*0.43);

apmap.prerun=Map(velposap.prerun(:,1:2),velposap.prerun(:,[1 3]),'smooth',1);
apmap.run=Map(velposap.run(:,1:2),velposap.run(:,[1 3]),'smooth',1);
safemap.prerun=Map(velpossafe.prerun(:,1:2),velpossafe.prerun(:,[1 3]),'smooth',1);
safemap.run=Map(velpossafe.run(:,1:2),velpossafe.run(:,[1 3]),'smooth',1);


safemap.postrun.all=Map(velpossafe.postrun.all(:,1:2),velpossafe.postrun.all(:,[1 3]),'smooth',1);
apmap.postrun.all=Map(velposap.postrun.all(:,1:2),velposap.postrun.all(:,[1 3]),'smooth',1);
safemap.postrun.firsts=Map(velpossafe.postrun.firsts(:,1:2),velpossafe.postrun.firsts(:,[1 3]),'smooth',1);
apmap.postrun.firsts=Map(velposap.postrun.firsts(:,1:2),velposap.postrun.firsts(:,[1 3]),'smooth',1);

figure('Position',[640 309 1167 626]);
subplot(3,2,1)
        plot(apmap.prerun.x,apmap.prerun.z,'r');
        Ylim=get(gca,'YLim');
%          if strcmp(prev.airpuff.dir,'LtoR')
%              arrow([prevairP-0.1 Ylim(2)*0.8],[prevairP Ylim(2)*0.8],'Length',5);
%          elseif strcmp(prev.airpuff.dir,'RtoL')
%              arrow([prevairP+0.1 Ylim(2)*0.8],[prevairP Ylim(2)*0.8],'Length',5);
%          end
        ylim([0 100]);
        PlotHVLines(airP,'v','r');
%          PlotHVLines(prevairP,'v','g');
        xlabel('airpuff laps - prerun' )
subplot(3,2,3);hold on;
        plot(apmap.run.x,apmap.run.z,'r');
        if strcmp(airpuff.dir,'LtoR')
            arrow([airP-0.1 Ylim(2)*0.8],[airP Ylim(2)*0.8],'Length',5);
        elseif strcmp(airpuff.dir,'RtoL')
            arrow([airP+0.1 Ylim(2)*0.8],[airP Ylim(2)*0.8],'Length',5);
        end
        ylim([0 100]);
        PlotHVLines(airP,'v','r');
        xlabel('airpuff laps - RUN')
subplot(3,2,5);hold on;
        plot(apmap.postrun.all.x,apmap.postrun.all.z,'k');
        plot(apmap.postrun.firsts.x,apmap.postrun.firsts.z,'r');
        ylim([0 100]);
        PlotHVLines(airP,'v','r');
        xlabel('airpuff laps - postrun');
subplot(3,2,2)
        plot(safemap.prerun.x,safemap.prerun.z);
%          if strcmp(prev.airpuff.dir,'LtoR')
%              arrow([prevairP-0.1 Ylim(2)*0.8],[prevairP Ylim(2)*0.8],'Length',5);
%          elseif strcmp(prev.airpuff.dir,'RtoL')
%              arrow([prevairP+0.1 Ylim(2)*0.8],[prevairP Ylim(2)*0.8],'Length',5);
%          end
        ylim([0 100]);
%          PlotHVLines(prevairP,'v','g');
        PlotHVLines(airP,'v','r');
        xlabel('safe laps - prerun');
subplot(3,2,4)
        ylim([0 100]);
        plot(safemap.run.x,safemap.run.z);
        PlotHVLines(airP,'v','r');
        xlabel('safe laps - run');
subplot(3,2,6);hold on;
        ylim([0 100]);
        plot(safemap.postrun.all.x,safemap.postrun.all.z,'k');
        plot(safemap.postrun.firsts.x,safemap.postrun.firsts.z,'b');
        PlotHVLines(airP,'v','r');
        xlabel('safe laps - postrun');
suptitle([xml '-AvergageVelocitites on laps'])

%  save([xml '-VelocityCurves.mat'],'safemap','apmap','airP','prevairP');
save([xml '-VelocityCurves.mat'],'safemap','apmap','airP');
cd('/media/Data-01/All-Rats/AllRats-VelocityCurves/')
plot2svg([xml '-VelocityCurves.svg'],gcf);


%  		    PlotColorIntervals([airP-0.01 airP+0.01],'Red','v');
%  		    xlabel(['RUN - Field peak: ' int2str(stats.peak) 'Hz. Field mean rate: ' int2str(stats.mean) 'Hz']);
%  		    % Plot Airpuff direction arrow
%  		    ylim=get(gca,'YLim');
		   