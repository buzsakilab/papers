function [LtoRlaps, RtoLlaps, Uturnlaps, leftLimit, rightLimit] = ClassifyLaps(varargin)

%ClassifyLaps - Classify laps as Left to Right, Right to Left and Uturnlaps.
%
%  USAGE
%
%    [LtoRlaps, RtoLlaps, Uturnlaps] = ClassifyLaps(runpos)
%
%    runs           run subsession numbers
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'session'		session path
%     'save'     	'on','off' (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    LtoRlaps       list of time intervals for Left to Right laps
%    RtoLlaps       list of time intervals for right to left laps
%    Uturnlaps      list of time intervals for U-turns
%
%  NOTE
%
%    lap distance set up to 150 cm (center of maze excluding ends)
%
%  SEE
%
%   CheckLaps.m
%
% Gabrielle Girardeau, December 2014
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
		case 'save',
			savevar= lower(varargin{i+1});
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

if ~isempty(session)
  cd (session)
  xml=session(end-13:end)
%    SetCurrentSession('filename',xml,'spikes','off')
  SetCurrentSession(xml,'spikes','off')
end

load('airpuff.mat');
pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);

load('runintervals.mat')
runint=runintervals;

runpos=Restrict(pos,runint);

LtoRlaps=[];
RtoLlaps=[];
Uturnlaps=[];

figure;
a=subplot(2,1,1);
%  PlotXY(runpos(:,2:3),'.');
PlotXY(runpos(:,2:3));
PlotHVLines(airpuff.loc*0.43,'v','r');
set(gca,'YDir','Reverse');
b=subplot(2,1,2);
PlotXY(runpos(:,1:2));
set(gcf,'Position',[355 106 1370 871]);
PlotHVLines(airpuff.loc*0.43,'h','r');
input('Define Left Limit')
Lpoint=get(gca,'CurrentPoint');
leftLimit=round(Lpoint(1,1));

%  input('Define Right Limit')
%  Rpoint=get(gca,'CurrentPoint');
%  rightLimit=round(Rpoint(1,1))
rightLimit=leftLimit+150;

lapindex=runpos(:,2)>leftLimit&runpos(:,2)<rightLimit;
lapIntervals=ToIntervals(runpos(:,1),lapindex);

lapIntervals(lapIntervals(:,2)-lapIntervals(:,1)>7000,:)=[];

for i=1:size(lapIntervals,1)
  currentInterval=lapIntervals(i,:);
  lappos=Restrict(runpos,currentInterval);
  lapbeg=lappos(1,2);
  lapend=lappos(end,2);
  if abs(lapbeg-lapend)>10
    if lapbeg>lapend
      RtoLlaps=[RtoLlaps;currentInterval];
    elseif lapbeg<lapend
      LtoRlaps=[LtoRlaps;currentInterval];
    end
  else
    if sum(abs(diff(lappos(:,2))))>80
      Uturnlaps=[Uturnlaps;lapIntervals(i,:)];
    end
  end
end

RtoLpos=Restrict(runpos,RtoLlaps);
LtoRpos=Restrict(runpos,LtoRlaps);
if ~isempty(Uturnlaps)
  Uturnpos=Restrict(runpos,Uturnlaps);
else
  Uturnpos=[];
end
subplot(b)
hold on

PlotXY(RtoLpos(:,1:2),'.g');
PlotXY(LtoRpos(:,1:2),'.r');
if ~isempty(Uturnpos)
  PlotXY(Uturnpos(:,1:2),'.k');
end
PlotColorIntervals(RtoLlaps,'PaleGreen','v');
PlotColorIntervals(LtoRlaps,'Coral','v');


if strcmp(savevar,'on') & ~isempty(session)
  save([xml '-Laps'],'RtoLlaps','LtoRlaps','Uturnlaps','rightLimit','leftLimit');
elseif strcmp(savevar,'on') & isempty(session)
  savename = input('Specify complete savename :','s')
  save(savename,'RtoLlaps','LtoRlaps','Uturnlaps','rightLimit','leftLimit');
end

