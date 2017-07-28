function [PrePmeanV,YPrePmeanV,ratio] = BehavioralMeasures (varargin)

%  BehavioralMeasures - Calculates Mean speed in the current/previous danger zone for preRUN, RUN and postRUN epochs and associated Ratios. Per session.
%
%  USAGE
%
%   [PrePmeanV,YPrePmeanV,ratio] = BehavioralMeasures (<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'session'		session path
%     'savevar'		save variables 'on'/'off' (default)
%    =========================================================================
%
%  OUTPUT
%
%   BehavioralMeasure per session (variable)
%
%   PrePmeanV   mean velocity in the danger zone before the current airpuff (.prerun, .run, .postrun)
%   YPrePmeanV  mean velocity in the danger zone before the previous airpuff (location of the previous training day).
%   Ratio       Velocity ratios of the current/previous dangerr zones (.prerun, .run, .postrun)
%
%  NOTE
%
%  SEE
%  
%    BehavioralMeasures_Stats, BehavioralMeasures_Plot
%
% Gabrielle Girardeau, March 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
session=[];
savevar='off';

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
		case 'session',
			session = varargin{i+1};
		case 'savevar',
			savevar = varargin{i+1};
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

session

if ~isempty(session)
  cd (session)
  xml=session(end-13:end)
  SetCurrentSession(xml,'spikes','off')
end

load('airpuff.mat');
if exist('airpuffPrevious.mat')==2
  Y=load('airpuffPrevious.mat');
else
  Y=[];
end

pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);

if exist('runintervals.mat')==2
  load('runintervals.mat');
  runint=runintervals;
end

%  runint=RunIntervals(runs);
prerun=runint(1,:);
run=runint(2,:);
postrun=runint(3,:);

lapfile=dir('*-Laps.mat');
load(lapfile.name);

allrunpos=Restrict(pos,runint);
vel=LinearVelocity(pos,5);
airpuffloc=airpuff.loc*0.43;
if ~isempty(Y)
  Yairpuffloc=Y.airpuff.loc*0.43;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% DAY AIRPUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%
zoneLeft=[airpuffloc-20 airpuffloc];
zoneRight=[airpuffloc airpuffloc+20];

Leftindex=allrunpos(:,2)>=zoneLeft(1)&allrunpos(:,2)<=zoneLeft(2);
Rightindex=allrunpos(:,2)>=zoneRight(1)&allrunpos(:,2)<zoneRight(2);

Leftintervals=ToIntervals(allrunpos(:,1),Leftindex);
Rightintervals=ToIntervals(allrunpos(:,1),Rightindex);

% Defining time intervals for PrePuff zone, PostPuff zone, PreSafe zone, PostSafe Zone.
if strcmp(airpuff.dir,'LtoR');
  statusPreP=InIntervals(Leftintervals(:,1),LtoRlaps);
  PrePint=Leftintervals(statusPreP,:);

  statusPostP=InIntervals(Rightintervals(:,1),LtoRlaps);
  PostPint=Rightintervals(statusPostP,:);

  statusPostS=InIntervals(Leftintervals(:,1),RtoLlaps);
  PostSint=Leftintervals(statusPostS,:);

  statusPreS=InIntervals(Rightintervals(:,1),RtoLlaps);
  PreSint=Rightintervals(statusPreS,:);
elseif strcmp(airpuff.dir,'RtoL');
  statusPreP=InIntervals(Rightintervals(:,1),RtoLlaps);
  PrePint=Rightintervals(statusPreP,:);

  statusPostP=InIntervals(Leftintervals(:,1),RtoLlaps);
  PostPint=Leftintervals(statusPostP,:);

  statusPostS=InIntervals(Rightintervals(:,1),LtoRlaps);
  PostSint=Rightintervals(statusPostS,:);

  statusPreS=InIntervals(Leftintervals(:,1),LtoRlaps);
  PreSint=Leftintervals(statusPreS,:);
end

% Subselecting per session
%%%%%% PreRun
prerunPrePint=ExcludeIntervals(PrePint,runint(2:3,:));
prerunPostPint=ExcludeIntervals(PostPint,runint(2:3,:));
prerunPreSint=ExcludeIntervals(PreSint,runint(2:3,:));
prerunPostSint=ExcludeIntervals(PostSint,runint(2:3,:));
prerunPrePvel=Restrict(vel,prerunPrePint);
PrePmeanV.prerun=nanmean(prerunPrePvel(:,2));

%%%%%% Run
runPrePint=ExcludeIntervals(PrePint,[runint(1,:);runint(3,:)]);
runPostPint=ExcludeIntervals(PostPint,[runint(1,:);runint(3,:)]);
runPreSint=ExcludeIntervals(PreSint,[runint(1,:);runint(3,:)]);
runPostSint=ExcludeIntervals(PostSint,[runint(1,:);runint(3,:)]);
if ~isempty(runPrePint)
  runPrePvel=Restrict(vel,runPrePint);
  PrePmeanV.run=nanmean(runPrePvel(:,2));
else
  PrePmeanV.run=NaN;
end

%%%%%% PostRun
postrunPrePint=ExcludeIntervals(PrePint,runint(1:2,:));
postrunPostPint=ExcludeIntervals(PostPint,runint(1:2,:));
postrunPreSint=ExcludeIntervals(PreSint,runint(1:2,:));
postrunPostSint=ExcludeIntervals(PostSint,runint(1:2,:));
if ~isempty(postrunPrePint)
  postrunPrePvel=Restrict(vel,postrunPrePint);
  PrePmeanV.postrun=nanmean(postrunPrePvel(:,2));
else
  PrePmeanV.postrun=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PREVIOUS DAY AIRPUFF %%%%%%%%%%%%%%%%%%%
if ~isempty(Y)
  YzoneLeft=[Yairpuffloc-20 Yairpuffloc];
  YzoneRight=[Yairpuffloc Yairpuffloc+20];

  YLeftindex=allrunpos(:,2)>=YzoneLeft(1)&allrunpos(:,2)<=YzoneLeft(2);
  YRightindex=allrunpos(:,2)>=YzoneRight(1)&allrunpos(:,2)<YzoneRight(2);

  YLeftintervals=ToIntervals(allrunpos(:,1),YLeftindex);
  YRightintervals=ToIntervals(allrunpos(:,1),YRightindex);

  % Defining time intervals for PrePuff zone, PostPuff zone, PreSafe zone, PostSafe Zone.
  Y.airpuff.dir
  Yairpuffloc
  airpuff.dir
  airpuffloc
  if strcmp(Y.airpuff.dir,'LtoR');
    YstatusPreP=InIntervals(YLeftintervals(:,1),LtoRlaps);
    YPrePint=YLeftintervals(YstatusPreP,:);
  elseif strcmp(Y.airpuff.dir,'RtoL');
    YstatusPreP=InIntervals(YRightintervals(:,1),RtoLlaps);
    YPrePint=YRightintervals(YstatusPreP,:);
  end

  % Subselecting per session
  %%%%%% PreRun
  YprerunPrePint=ExcludeIntervals(YPrePint,runint(2:3,:));

  YprerunPrePvel=Restrict(vel,YprerunPrePint);
  YPrePmeanV.prerun=nanmean(YprerunPrePvel(:,2));

  %%%%%% Run
  YrunPrePint=ExcludeIntervals(YPrePint,[runint(1,:);runint(3,:)]);
  if ~isempty(YrunPrePint)
    YrunPrePvel=Restrict(vel,YrunPrePint);
    YPrePmeanV.run=nanmean(YrunPrePvel(:,2));
  else
    YPrePmeanV.run=NaN;
  end
  %%%%%% PostRun
  YpostrunPrePint=ExcludeIntervals(YPrePint,runint(1:2,:));

  YpostrunPrePvel=Restrict(vel,YpostrunPrePint);
  YPrePmeanV.postrun=nanmean(YpostrunPrePvel(:,2));
else
  YPrePmeanV.prerun=NaN;
  YPrePmeanV.run=NaN;
  YPrePmeanV.postrun=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% RATIOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio.prerun=(YPrePmeanV.prerun-PrePmeanV.prerun)/(YPrePmeanV.prerun+PrePmeanV.prerun);
ratio.run=(YPrePmeanV.run-PrePmeanV.run)/(YPrePmeanV.run+PrePmeanV.run);
ratio.postrun=(YPrePmeanV.postrun-PrePmeanV.postrun)/(YPrePmeanV.postrun+PrePmeanV.postrun);

%%%%% OVERLAP ?
if ~ isempty(Y)
  if strcmp(Y.airpuff.dir,airpuff.dir) & abs(Yairpuffloc-airpuffloc)<20
    overlap=1;
  else
    overlap=0;
  end
else
  overlap=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1:2)
%  PlotXY(vel,'color',rgb('ForestGreen'));
PlotXY(vel(1:15:end,:),'color',rgb('ForestGreen'));
hold on;
%  PlotXY(pos(:,1:2),'color',rgb('DimGray'));
PlotXY(pos(1:15:end,1:2),'color',rgb('DimGray'));
PrePpos=Restrict(pos,PrePint);
PlotXY(PrePpos,'.','color',rgb('Crimson'));
PlotHVLines(airpuffloc,'h','color',rgb('Crimson'));
if ~isempty(Y)
  YPrePpos=Restrict(pos,YPrePint);
  PlotXY(YPrePpos,'.','color',rgb('MidnightBlue'));
  PlotHVLines(Yairpuffloc,'h','color',rgb('MidnightBlue'));
  PlotColorIntervals([YzoneLeft;YzoneRight;zoneLeft;zoneRight],'LightGray','h');
else
  PlotColorIntervals([zoneLeft;zoneRight],'LightGray','h');
end

ylim([0 300]);
if ~isempty(session)
  xlabel(xml);
end

subplot(2,2,3)
bar([YPrePmeanV.prerun PrePmeanV.prerun 0 YPrePmeanV.run PrePmeanV.run 0 YPrePmeanV.postrun PrePmeanV.postrun],'FaceColor',rgb('MediumSeaGreen'),'LineStyle','none');
set(gca,'xTickLabel',{'Y-prerun' 'Day-prerun' '' 'Y-run' 'run' '' 'Y-postrun' 'day-postrun'});

if overlap
  ylabel('!! OVERLAPPING AIRPUFF LOC/DIR')
end

subplot(2,2,4)
bar([ratio.prerun ratio.run ratio.postrun],'FaceColor',rgb('SeaGreen'),'LineStyle','none');
set(gca,'xTickLabel',{'Ratio-prerun' 'Ratio-run' 'Ratio-postrun'});
xlabel('Ratios');

set(gcf,'Position',[546 73 1194 889]);

plot2svg(['/media/Data-01/All-Rats/BehavioralMeasures/' xml '-BehavioralMeasures.svg'],gcf);
print(gcf,['/media/Data-01/All-Rats/BehavioralMeasures/' xml '-BehavioralMeasures'],'-dpng');

if strcmp(savevar,'on')
  save('BehavioralMeasures.mat','PrePmeanV','YPrePmeanV','ratio');
end
