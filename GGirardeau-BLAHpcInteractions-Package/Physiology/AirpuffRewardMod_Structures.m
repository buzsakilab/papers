function [AllCurves,CellIndex,x,FRIndex,AllAPpeth,AllRWpeth] = AirpuffRewardMod_Structures(structure,varargin)


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

load('/media/Data-01/All-Rats/sessionindexing.mat')
load('/media/Data-01/All-Rats/Structures/structures.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/SpikeParameters.mat');

struct=eval(structure);

AllCurves.all=[];
AllCurves.airpuff=[];
AllCurves.safe=[];
AllCurves.prerun=[];
AllCurves.postrun=[];


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
	airpuffcurve.safe=flip(MazeCurve{cells(j)}.run.safelapscurve.rate(APzone));
	airpuffcurve.prerun=MazeCurve{cells(j)}.prerun.curve.rate(APzone);
	airpuffcurve.postrun=MazeCurve{cells(j)}.postrun.curve.rate(APzone);
      elseif strcmp(Airpuff.dir,'RtoL')
	airpuffcurve.all=flip(MazeCurve{cells(j)}.run.curve.rate(APzone));
	airpuffcurve.airpuff=flip(MazeCurve{cells(j)}.run.airpufflapscurve.rate(APzone));
	airpuffcurve.safe=(MazeCurve{cells(j)}.run.safelapscurve.rate(APzone));
	airpuffcurve.prerun=MazeCurve{cells(j)}.prerun.curve.rate(APzone);
	airpuffcurve.postrun=MazeCurve{cells(j)}.postrun.curve.rate(APzone);
      end
      
      
      if sum(airpuffcurve.all)~=0 & sum(airpuffcurve.airpuff)~=0 & sum(airpuffcurve.safe)~=0
	AllCurves.all=[AllCurves.all;airpuffcurve.all];
	AllCurves.airpuff=[AllCurves.airpuff;airpuffcurve.airpuff];
	AllCurves.safe=[AllCurves.safe;airpuffcurve.safe];
	AllCurves.prerun=[AllCurves.prerun;airpuffcurve.prerun];
	AllCurves.postrun=[AllCurves.postrun;airpuffcurve.postrun];
	CellIndex=[CellIndex;ratsess structurecells(structurecells(:,3)==cells(j),:)];
	FRIndex=[FRIndex;cellFR(j)];
	
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
  save(['/media/Data-01/All-Rats/AllRats-AirpuffRewardColorPlots/' structure '.mat'],'AllCurves','CellIndex','x','FRIndex','AllAPpeth','tAP','AllRWpeth','tRW');
end


