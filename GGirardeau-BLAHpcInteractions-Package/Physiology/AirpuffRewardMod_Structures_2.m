function [AllCurves,CompleteCurves,CellIndex,x,FRIndex,AllAPpeth,AllRWpeth] = AirpuffRewardMod_Structures_2(structure,varargin)


% Defaults
ctype = 'pyr';
savevar = 'off';

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
    case 'ctype',
      ctype = varargin{i+1};
    case 'savevar'
      savevar = varargin{i+1};
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

load('/media/Data-01/All-Rats/sessionindexing.mat');
load('/media/Data-01/All-Rats/Structures/structures.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/SpikeParameters.mat');

struct=eval(structure);

AllCurves.all=[];
AllCurves.airpuff=[];
AllCurves.safe.flip=[];
AllCurves.safe.unflip=[];
AllCurves.prerun=[];
AllCurves.postrun=[];

CompleteCurves.run.all=[];
CompleteCurves.run.laps=[];
CompleteCurves.prerun.all=[];
CompleteCurves.postrun.all=[];
CompleteCurves.prerun.laps=[];
CompleteCurves.postrun.laps=[];
CompleteCurves.safe=[];
CompleteCurves.airpuff=[];

AllAPpeth.short=[];
AllAPpeth.long=[];

AllRWpeth.short=[];
AllRWpeth.long=[];

FRIndex=[];
CellIndex=[];


for i=1:length(xmlpath)
  xml=xmlpath{i}(end-14:end-1);
  cd(xmlpath{i});
  ratsess=ratsessionindex(strcmp(xmlpath{i},xmlpath),1:2); %[rat session]
  structureshanks=struct(ismember(struct(:,1:2),ratsess,'rows'),3);
  
  if exist([xml '-AirpuffRewardMod-SafeRew.mat'])==2 & exist([xml '-AirpuffRewardMod-2.mat'])==2 & ~isempty(structureshanks) & ~strcmp(xml,'Rat08-20130708')
    load([xml '-AirpuffRewardMod-2.mat']);
    load([xml '-AirpuffRewardMod-SafeRew.mat']);
                     
    sessiontypes=finalType(ismember(finalType(:,1:2),ratsess,'rows'),3:5);	% shank/cell/type for all session cells
    sessionFR=SpikeParameters(ismember(SpikeParameters(:,1:2),ratsess,'rows'),[3 4 10]);
   
    structurecellsID=Index(ismember(Index(:,1),structureshanks,'rows'),3);
    structureFR=sessionFR(ismember(sessionFR(:,1),structureshanks,'rows'),3);
    structurecells=Index(ismember(Index(:,1),structureshanks,'rows'),:);
    structuretypes=sessiontypes(ismember(sessiontypes(:,1),structureshanks,'rows'),3);
    
    if strcmp(ctype,'all')
      cells=structurecellsID;
      cellFR=structureFR;
    elseif strcmp(ctype,'pyr')
      cells=structurecellsID(structuretypes==1);
      cellFR=structureFR(structuretypes==1);
    end
           
    % Define APzone 10bins=20cm
    [~,airPbin]=find(abs(MazeCurve{cells(1)}.run.curve.x-Airpuff.loc)==min(abs(MazeCurve{cells(1)}.run.curve.x-Airpuff.loc)));
    APzone=zeros(1,length(MazeCurve{cells(1)}.run.curve.x));
    APzone(airPbin-20:airPbin+20)=1;
    APzone=logical(APzone);
    for j=1:length(cells)
      if strcmp(Airpuff.dir,'LtoR')
	airpuffcurve.all=MazeCurve{cells(j)}.run.curve.rate(APzone);
	airpuffcurve.airpuff=MazeCurve{cells(j)}.run.airpufflapscurve.rate(APzone);
	airpuffcurve.safe.flip=flip(MazeCurve{cells(j)}.run.safelapscurve.rate(APzone));	% safe flip is reverse to match approach to airpuff/airpuff location
	airpuffcurve.safe.unflip=MazeCurve{cells(j)}.run.safelapscurve.rate(APzone);		% safe unflip is maintained as absolute position so that pre-airpuff matches post-airpuff loc on safe laps
	airpuffcurve.prerun=MazeCurve{cells(j)}.prerun.curve.rate(APzone);
	airpuffcurve.postrun=MazeCurve{cells(j)}.postrun.curve.rate(APzone);
      elseif strcmp(Airpuff.dir,'RtoL')
	airpuffcurve.all=flip(MazeCurve{cells(j)}.run.curve.rate(APzone));
	airpuffcurve.airpuff=flip(MazeCurve{cells(j)}.run.airpufflapscurve.rate(APzone));
	airpuffcurve.safe.unflip=flip(MazeCurve{cells(j)}.run.safelapscurve.rate(APzone));
	airpuffcurve.safe.flip=MazeCurve{cells(j)}.run.safelapscurve.rate(APzone);
	airpuffcurve.prerun=MazeCurve{cells(j)}.prerun.curve.rate(APzone);
	airpuffcurve.postrun=MazeCurve{cells(j)}.postrun.curve.rate(APzone);
      end
      curve.prerun.all=MazeCurve{cells(j)}.prerun.curve.rate;
      curve.prerun.laps=MazeCurve{cells(j)}.prerunlaps.curve.rate(MazeCurve{cells(j)}.prerunlaps.curve.time~=0);
      curve.postrun.all=MazeCurve{cells(j)}.postrun.curve.rate;
      curve.postrun.laps=MazeCurve{cells(j)}.postrunlaps.curve.rate(MazeCurve{cells(j)}.postrunlaps.curve.time~=0);
      curve.run.all=MazeCurve{cells(j)}.run.curve.rate;
      
      curve.run.laps=MazeCurve{cells(j)}.run.alllapscurve.rate(MazeCurve{cells(j)}.run.alllapscurve.time~=0);
      curve.run.airpuff=MazeCurve{cells(j)}.run.airpufflapscurve.rate(MazeCurve{cells(j)}.run.airpufflapscurve.time~=0);
      curve.run.safe=MazeCurve{cells(j)}.run.safelapscurve.rate(MazeCurve{cells(j)}.run.safelapscurve.time~=0);
      
      if length(curve.prerun.laps)>76 | length(curve.postrun.laps)>76
	warning('Curve length fuckup')
      end
      curve.prerun.laps=curve.prerun.laps(1:75);
      curve.postrun.laps=curve.postrun.laps(1:75);
            
      if length(curve.run.airpuff)>76
	xml
	length(curve.run.airpuff)
	length(curve.run.safe)
      end
      if length(curve.run.safe)==76 | length(curve.run.airpuff)==76
	curve.run.safe=curve.run.safe(1:75);
	curve.run.laps=curve.run.laps(1:75);
	curve.run.airpuff=curve.run.airpuff(1:75);
      end
      
%        if sum(airpuffcurve.all)~=0 & sum(airpuffcurve.airpuff)~=0 & sum(airpuffcurve.safe.flip)~=0
      if nSpikes{cells(j)}.runalllaps>50 | nSpikes{cells(j)}.prerun.laps>50
	AllCurves.all=[AllCurves.all;airpuffcurve.all];
	AllCurves.airpuff=[AllCurves.airpuff;airpuffcurve.airpuff];
	AllCurves.safe.flip=[AllCurves.safe.flip;airpuffcurve.safe.flip];
	AllCurves.safe.unflip=[AllCurves.safe.unflip;airpuffcurve.safe.unflip];
	AllCurves.prerun=[AllCurves.prerun;airpuffcurve.prerun];
	AllCurves.postrun=[AllCurves.postrun;airpuffcurve.postrun];
	CellIndex=[CellIndex;ratsess structurecells(structurecells(:,3)==cells(j),:)];
	FRIndex=[FRIndex;cellFR(j)];
	
	CompleteCurves.run.all=[CompleteCurves.run.all;curve.run.all];
	CompleteCurves.run.laps=[CompleteCurves.run.laps;curve.run.laps];
	CompleteCurves.prerun.all=[CompleteCurves.prerun.all;curve.prerun.all];
	CompleteCurves.prerun.laps=[CompleteCurves.prerun.laps;curve.prerun.laps];
	CompleteCurves.postrun.all=[CompleteCurves.postrun.all;curve.postrun.all];
	CompleteCurves.postrun.laps=[CompleteCurves.postrun.laps;curve.postrun.laps];
	CompleteCurves.safe=[CompleteCurves.safe;curve.run.safe];
	CompleteCurves.airpuff=[CompleteCurves.airpuff;curve.run.airpuff];
	
	[airpuff.short,t.short]=SyncHist(AirpuffPETH{cells(j)}.raster,AirpuffPETH{cells(j)}.indices,'mode','sum','smooth',0);
	if isnan(airpuff.short)
	  airpuff.short=zeros(100,1);
	end
	[airpuff.long,t.long]=SyncHist(AirpuffPETH{cells(j)}.raster,AirpuffPETH{cells(j)}.indices,'mode','sum','smooth',0,'durations',[-2 2]);
	if isnan(airpuff.long)
	  airpuff.long=zeros(100,1);
	end
	if isempty(airpuff.long)
	  airpuff.short=zeros(100,1);
	  airpuff.long=zeros(100,1);
	end
	AllAPpeth.short=[AllAPpeth.short;airpuff.short'];
	AllAPpeth.long=[AllAPpeth.long;airpuff.long'];
	
	[reward.short,t.short]=SyncHist(SafeRewPETH{cells(j)}.run.raster,SafeRewPETH{cells(j)}.run.indices,'mode','sum','smooth',0);
	if isnan(reward.short)
	  reward.short=zeros(100,1);
	end
	[reward.long,t.long]=SyncHist(SafeRewPETH{cells(j)}.run.raster,SafeRewPETH{cells(j)}.run.indices,'mode','sum','smooth',0,'durations',[-2 2]);
	if isnan(reward.long)
	  reward.long=zeros(100,1);
	end
	AllRWpeth.short=[AllRWpeth.short;reward.short'];
	AllRWpeth.long=[AllRWpeth.long;reward.long'];
      end
    end
    tAP.short=t.short;
    tAP.long=t.long;
    tRW.short=t.short;
    tRW.long=t.long;  
  end
end
x=-40:2:40;

if strcmp(savevar,'on')
  save(['/media/Data-01/All-Rats/AllRats-AirpuffRewardColorPlots/' structure '-' ctype '.mat'],'AllCurves','CompleteCurves','CellIndex','x','FRIndex','AllAPpeth','tAP','AllRWpeth','tRW');
end


