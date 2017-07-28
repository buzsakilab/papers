function [FinalExcMonoSynList,FinalInhMonoSynList,ccgR,tR,Pred,FinalExcCellList,FinalInhCellList,oops,completeIndex,Bounds] = MonoSynConvClick (binSize,duration,alpha,varargin)

%MonoSynConvClick - Finds all monosynaptic connections for a given recording day (all shanks, all cells/clu) using the convolution method (Eran Stark), with manual, clickable additional sorting.
%
%  USAGE
%
%    [ExcMonoS,InhMonoS,ccgR,ExcCellList,InhCellList] = MonoSyn (binSize,duration,alpha,<options>)
%
%    binSize            size of the CCG bin
%    duration           total duration of CCG in ms (if 40 : from -20 to 20)
%    alpha		significance level (0.001,0.01 or 0.05)
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'manual'		Enable manual suppression of spurious Excitatory/Inhibitory
%  			cells and monosynaptic connections
%  			'on' or 'off' (default = 'off')
%     'show'     	Show final Monosynaptic connection for defined cell types
%  			'on' or 'off' (default = 'off')
%     'session'		session : e.g. 'Rat08-20130717'
%     'DB'		DataBase name if figures/variables are to be stored in DB (default = [])
%     'folder'		Name of folder where to save the png figures. Default 'none'. Folder must be pre-created in session folder.
%    =========================================================================
%
%  OUTPUT
%
%    ExcMonoS       	List of all excitatory monosynaptic connections : shank(cell1) / cell# / shank(cell2) / cell#
%    InhMonoS		Same for inhibitory connections
%    ccgR           	CCG 3D matrix for all cells (nbins x ncells x ncells)
%    ExcCellList	List of all excitatory cells : shank / cell#
%    InhCellList	List of all inhibitory cells : shank / cell#
%    oops		List of cells that were found to be both inhibitory and excitatory (usually bad clusters with shared noise)
%    completeIndex	List of all cells with associated ID (to retrieve specific data in ccgR, for ex.) : shank/cell#/ID
%
%  EXAMPLE
%
%    [FinalExcM,FinalInhM,ccgR,tR,Pred,FinalExcCellList,FinalInhCellList,oops,Idx,Bounds]=MonoSynConvClick(0.5,80,0.001,'manual','on','show','on','session','Rat08-20130708','DB','Rat08_MonoSynapticConnections');
%
%  NOTES
%
%    !!! Session must be loaded BEFORE using MonoSynConvClick using FMAToolbox's SetCurrentSession !!!
%
%    The Bonferroni correction is applied to singifiance thresold (alpha) to compensate for multiple comparisons.
%    The effective alpha has to be entered as parameter (0.001, 0.01 or 0.05) but the significance tests will be performed with alpha/nBIns (hardcoded for now!)
%
%    TEST WINDOW (1ms<=window<4ms) AND BONFERRONI CORRECTION (n=6 bins) HARDCODED
%    CONVOLUTION WINDOW HARDCODED (10 ms)
%
%  SEE
%
%    See also : CCG.m (from FMAToolbox), cch_conv, AxesInsetBars (Brendon Watson) (dependencies)
%
% Jan 2014 by Gabrielle Girardeau calling Eran Stark's cch_conv.m, FMAToolbox CCG.m
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
show = 'off' ;
manual = 'off';
folder = 'none';

global session
global db;
db = [];
session = [];

% Check number of inputs
if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
end

% Check varargin
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters  (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'show',
			show = varargin{i+1};
			if ~isstring(show),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
			end
		case 'manual',
			manual = varargin{i+1};
			if ~isstring(manual),
				error('Incorrect value for property ''manual'' (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
			end
		case 'db',
			db = varargin{i+1};
			if ~isstring(db),
				error('Incorrect value for property ''db'' (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
			end
		case 'session',
			session = varargin{i+1};
			if ~isstring(session),
				error('Incorrect value for property ''db'' (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).');
			end
		case 'folder'
			folder = varargin{i+1};
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MonoSyn">MonoSyn</a>'' for details).']);
	end
end

nBins=duration/binSize+1;

if ~isempty(session)
  cd (session)
  xml=session(end-13:end)
  SetCurrentSession(xml)
else
  xml=input('Enter session basename - no brackets','s')
end

% Warning if variable already exists
if exist([xml '-MonoSynConvClick.mat'])==2;
  keyin = input('Variable with MonoSynConvClick output already exists. Running the function will overwrite it : do you wish to proceed? y/n','s');
  keypress=1;
  while keypress
    if strcmp(keyin,'n')
      error('Aborted')
    elseif strcmp(keyin,'y')
      keypress=0;
    end
  end
end

if strcmp(folder,'none')==0
  if exist(folder)~=7
    mkdir(folder)
  end
end

%  Create allspikes with : spiketimes / shank / cluster(cell#) / uniqueID
%  MUA and artefacts are suppressed wich means : cell# 0 and 1 do not exist --> renumber IDs from 1 to n as needed by CCG.m
spikes1=GetSpikeTimes('output','full');
spikes2=GetSpikeTimes('output','numbered');
allspikes=[spikes1 spikes2(:,2)];
% Remove noise and MUA clusters
allspikes(allspikes(:,3)==0,:)=[];
allspikes(allspikes(:,3)==1,:)=[];

%  %  Restrict cells for function testing
%  allspikes(allspikes(:,4)>20,:)=[];
%  warning('Number of cells restricted for testing');

%  Renumber IDs from 1 to n as needed by CCG.m
unik=unique(allspikes(:,4));
for ii=1:length(unik)
    allspikes(allspikes(:,4)==unik(ii),4)=ii;
end
if sum(diff(unique(allspikes(:,4))))~=length(unik)-1
    disp('ID renumbering failed')
    return
end

% convert binSize and duration to seconds (spiketimes are in s)
binSize=binSize/1000;
duration=duration/1000;

spikesIDs = allspikes(:,4); 	% spikes unique IDs only
spiketimes = allspikes(:,1);	% spikes times only
IDindex = unique(spikesIDs);	% list of IDs

% Make a complete index of all cells : shank/cell/ID : te be used later to identify cells from ID and vice-versa
completeIndex=unique(allspikes(:,2:4),'rows');
completeIndex=sortrows(completeIndex,3);

% Initialize lists
global SpuriousExcMonoSynList;
global SpuriousInhMonoSynList;
SpuriousExcMonoSynList=[];
SpuriousInhMonoSynList=[];


% Create CCGs (including autoCG) for all cells
[ccgR,tR] = CCG(spiketimes,spikesIDs,'binSize',binSize,'duration',duration);
% Initialize MonoSynaptic Lists and other 3D result matrices
ExcMonoS=[];
InhMonoS=[];
Pred=zeros(length(tR),length(IDindex),length(IDindex));
Bounds=zeros(2,length(IDindex),length(IDindex));

for refcellID=1:length(IDindex)
  for cell2ID=1:length(IDindex)
    if refcellID~=cell2ID 				% exclude autocorrelos
      cch=ccgR(:,refcellID,cell2ID);			% extract corresponding cross-correlation histogram vector
      refcellshank=completeIndex(completeIndex(:,3)==refcellID,1);
      cell2shank=completeIndex(completeIndex(:,3)==cell2ID,1);
      if refcellshank==cell2shank
	% Calculate predictions using Eran's cch_conv correcting for zero-lag bins (cch_conv without zerolag bin then replaced with NaN) on same-shank cells !! HARDCODED for BinSize=0.5ms (3 central bins ignored)
	sameshankcch=cch;
	sameshankcch(floor(nBins/2):ceil(nBins/2)+1)=[];
	[pvals,pred,qvals]=cch_conv(sameshankcch,20);
	pred=[pred(1:length(pred)/2);NaN;NaN;NaN;pred(length(pred)/2+1:end)];
      else
	% calculate predictions using Eran's cch_conv
	[pvals,pred,qvals]=cch_conv(cch,20);
      end
      % Store predicted values and pvalues for subsequent plotting
      Pred(:,refcellID,cell2ID)=pred;
      % window for monosyn connections
      window=tR>=0.001&tR<0.004;		%%%%%%%%%%% HARDCODED %%%%%%%%%%%%
      nBonf=6;					%%%%%%%%%%% HARDCODED %%%%%%%%%%%% (Bonferroni correction for test on 6 bins)
      % Calculate upper and lower limits
      hiBound=poissinv(1-alpha/nBonf,max(pred(window)));
      loBound=poissinv(alpha/nBonf, min(pred(window)));
      Bounds(1,refcellID,cell2ID)=hiBound;
      Bounds(2,refcellID,cell2ID)=loBound;
      % Find if significant periods falls in monosynaptic window - Make a list
      % of all significant excitatory and inhibitory connections (IDs)
      sigExc=cch>hiBound;
      sigInh=cch<loBound;
      if sum(sigExc(window))~=0
	ExcMonoS=[ExcMonoS;refcellID cell2ID];
      end
      if sum(sigInh(window))~=0
	if ~(sigInh(tR==0.001)==1&&sum(sigInh(window))==1)
	  InhMonoS=[InhMonoS;refcellID cell2ID];
	end
      end
    end
  end
end

%  if ~isempty(db);
%      DBConnect
%      DBUse(db);
%      DBAddVariable(ccgR,['Hpc-Amygdala'],[session ' CCG-AllCells CCG'],'','','MonoSynConv.m');
%      DBAddVariable(tR,['Hpc-Amygdala'],[session ' CCG-AllCells TimeBins'],'','','MonoSynConv.m');
%      DBAddVariable(Pred,['Hpc-Amygdala'],[session ' CCG-AllCells Pred'],'','','MonoSynConv.m');
%  end



% Make a list of all excitatory cells (IDs)
if ~isempty(ExcMonoS)
    ExcCellList=unique(ExcMonoS(:,1));
else
    ExcCellList=[];
end
% Make a list of all inhibitory cells (IDs)
if ~isempty(InhMonoS)
    InhCellList=unique(InhMonoS(:,1));
else
    InhCellList=[];
end


% Manual sorting
if strcmp(manual,'on')
  position = [574 175 1072 730];
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:size(ExcCellList,1)
      refcell=ExcCellList(i);
      excitedcells = ExcMonoS(ExcMonoS(:,1)==refcell,2);
      refshankcell=completeIndex(completeIndex(:,3)==refcell,1:2);
      fg = figure('Position',position,'CloseRequestFcn',@CloseFigureFcnExc);
      suptitle([xml ' - Putative Excitatory Cell : shank ' int2str(refshankcell(1)) ' cell ' int2str(refshankcell(2)) ' (ID : ' int2str(refcell) ')']);
      for j=1:size(excitedcells,1)
	acellinfo.cell = [refcell excitedcells(j)];
	acellinfo.type = 'exc';
	a = SquareSubplot(size(excitedcells,1)+1,j);
	SpuriousExcMonoSynList(end+1,:)=acellinfo.cell;
	bar(tR,(ccgR(:,refcell,excitedcells(j))));
	hold on;
	% Plot predicted values
	plot(tR,Pred(:,refcell,excitedcells(j)),'g');
	%Plot upper and lower boundaries
	PlotHVLines(Bounds(1,refcell,excitedcells(j)),'h','r--');
	PlotHVLines(Bounds(2,refcell,excitedcells(j)),'h','r--');
	% Plot signif increased bins in red
	exc=ccgR(:,refcell,excitedcells(j));
	exc(exc<Bounds(1,refcell,excitedcells(j))|~window)=0;
	bar(tR,exc,'r');
	% Plot signif lower bins in blue
	inh=ccgR(:,refcell,excitedcells(j));
	inh(inh>Bounds(2,refcell,excitedcells(j))|~window)=0;
	bar(tR,inh,'c');
	xlim([-0.02 0.02]);
	PlotColorIntervals([0.00075 0.00375],'LightGrey','v');
	excitedshankcell=completeIndex(completeIndex(:,3)==excitedcells(j),1:2);
	xlabel(['shank ' int2str(excitedshankcell(1)) ' cell ' int2str(excitedshankcell(2)) ' (ID:' int2str(excitedcells(j)) ')' ]);
	set(a,'UserData',acellinfo,'Color',[1 .75 .75],'ButtonDownFcn',@subplotclick);
	% Plot an inset with the excitedcell ACG
	thisacg = ccgR(:,excitedcells(j),excitedcells(j));
	axh = AxesInsetBars(a,.2,[.5 .5 .5],tR',thisacg);
	axhpos = get(axh,'Position');
	set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)],'XTickLabel',[]);
      end
      SquareSubplot(size(excitedcells,1)+1,size(excitedcells,1)+1);
      bar(tR,ccgR(:,refcell,refcell),'k');
      xlim([-0.04 0.04]);
      xlabel('Reference Cell ACG');
      guidata(fg,refshankcell);
      waitfor(fg)
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:size(InhCellList,1)
      refcell=InhCellList(i);
      inhibitedcells = InhMonoS(InhMonoS(:,1)==refcell,2);
      refshankcell=completeIndex(completeIndex(:,3)==refcell,1:2);
      fg = figure('Position',position,'CloseRequestFcn',@CloseFigureFcnInh);
      suptitle([xml ' - Putative Inhibitory Cell : shank ' int2str(refshankcell(1)) ' cell ' int2str(refshankcell(2)) ' (ID : ' int2str(refcell) ')']);
      for j=1:size(inhibitedcells,1)
      	acellinfo.cell = [refcell inhibitedcells(j)];
	acellinfo.type = 'inh';
	a = SquareSubplot(size(inhibitedcells,1)+1,j);
	SpuriousInhMonoSynList(end+1,:)=acellinfo.cell;
	bar(tR, ccgR(:,refcell,inhibitedcells(j)));
	hold on;
	plot(tR,Pred(:,refcell,inhibitedcells(j)),'g');
	%Plot upper and lower boundaries
	PlotHVLines(Bounds(1,refcell,inhibitedcells(j)),'h','r--');
	PlotHVLines(Bounds(2,refcell,inhibitedcells(j)),'h','r--');
	% Plot signif increased bins in red
	exc=ccgR(:,refcell,inhibitedcells(j));
	exc(exc<Bounds(1,refcell,inhibitedcells(j))|~window)=0;
	bar(tR,exc,'r');
	% Plot signif lower bins in blue
	inh=ccgR(:,refcell,inhibitedcells(j));
	inh(inh>Bounds(2,refcell,inhibitedcells(j))|~window)=0;
	bar(tR,inh,'c');
	xlim([-0.02 0.02]);
	PlotColorIntervals([0.00075 0.00375],'LightGrey','v');
	inhibshankcell=completeIndex(completeIndex(:,3)==inhibitedcells(j),1:2);
	xlabel(['shank ' int2str(inhibshankcell(1)) ' cell ' int2str(inhibshankcell(2)) ' (ID:' int2str(inhibitedcells(j)) ')' ]);
	set(a,'UserData',acellinfo,'Color',[1 .75 .75],'ButtonDownFcn',@subplotclick);
	% Plot an inset with the inhibitedcell ACG
	thisacg = ccgR(:,inhibitedcells(j),inhibitedcells(j));
	axh = AxesInsetBars(a,.2,[.5 .5 .5],tR',thisacg);
	axhpos = get(axh,'Position');
	set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)],'XTickLabel',[]);
      end
      SquareSubplot(size(inhibitedcells,1)+1,size(inhibitedcells,1)+1);
      bar(tR,ccgR(:,refcell,refcell),'k');
      xlim([-0.04 0.04]);
      xlabel('Reference Cell ACG');
      guidata(fg,refshankcell);
      waitfor(fg)
  end
end

%  if ~isempty(db)
%    DBConnect
%    DBUse(db);
%    DBAddVariable(ExcCellList,['Hpc-Amygdala'],[session ' ExcCellList'],'','','MonoSynConv.m');
%    DBAddVariable(InhCellList,['Hpc-Amygdala'],[session ' InhCellList'],'','','MonoSynConv.m');
%    DBAddVariable(oops,['Hpc-Amygdala'],[session ' Oops'],'','','MonoSynConv.m');
%    DBAddVariable(ExcMonoS,['Hpc-Amygdala'],[session ' ExcMonoS'],'','','MonoSynConv.m');
%    DBAddVariable(InhMonoS,['Hpc-Amygdala'],[session ' InhMonoS'],'','','MonoSynConv.m');
%    DBAddVariable(completeIndex,['Hpc-Amygdala'],[session ' Cell Index'],'','','MonoSynConv.m');
%  end

disp('Done reviewing the connections - Calculating final lists. Final results will be optionnally plotted and saved.');

% Update Final Exc and Inh Lists
if ~isempty(SpuriousExcMonoSynList)
  FinalExcMonoSynList=ExcMonoS(~ismember(ExcMonoS,SpuriousExcMonoSynList,'rows'),:);
else
  FinalExcMonoSynList=ExcMonoS;
end
if ~isempty(SpuriousInhMonoSynList)
  FinalInhMonoSynList=InhMonoS(~ismember(InhMonoS,SpuriousInhMonoSynList,'rows'),:);
else
  FinalInhMonoSynList=InhMonoS;
end
FinalExcCellList=unique(FinalExcMonoSynList(:,1));
FinalInhCellList=unique(FinalInhMonoSynList(:,1));
oops=intersect(FinalExcCellList,FinalInhCellList);
if ~isempty(oops)
  warning('Some cells are still classified as both excitatory and inhibitory')
end


%  %  Plotting final results
if strcmp(show,'on')
  for i=1:size(FinalExcCellList,1)
    refcell=FinalExcCellList(i);
    excitedcells = FinalExcMonoSynList(FinalExcMonoSynList(:,1)==refcell,2);
    refshankcell=completeIndex(completeIndex(:,3)==refcell,1:2);
    fg=figure;
    suptitle([xml ' - Putative Excitatory Cell : shank ' int2str(refshankcell(1)) ' cell ' int2str(refshankcell(2)) ' (ID : ' int2str(refcell) ')']);
    for j=1:size(excitedcells,1)
      a=SquareSubplot(size(excitedcells,1)+1,j);
      bar(tR,(ccgR(:,refcell,excitedcells(j))));
      hold on;
      % Plot predicted values
      plot(tR,Pred(:,refcell,excitedcells(j)),'g');
      % Plot upper and lower boundaries
      PlotHVLines(Bounds(1,refcell,excitedcells(j)),'h','r--');
      PlotHVLines(Bounds(2,refcell,excitedcells(j)),'h','r--');
      % Plot signif increased bins in red
      exc=ccgR(:,refcell,excitedcells(j));
      exc(exc<Bounds(1,refcell,excitedcells(j))|~window)=0;
      bar(tR,exc,'r');
      % Plot signif lower bins in blue
      inh=ccgR(:,refcell,excitedcells(j));
      inh(inh>Bounds(2,refcell,excitedcells(j))|~window)=0;
      bar(tR,inh,'c');
      xlim([-0.02 0.02]);
      PlotColorIntervals([0.00075 0.00375],'LightGrey','v');
      excshankcell=completeIndex(completeIndex(:,3)==excitedcells(j),1:2);
      xlabel(['shank ' int2str(excshankcell(1)) ' cell ' int2str(excshankcell(2)) ' (ID:' int2str(excitedcells(j)) ')' ]);
      % Plot an inset with the excited cell ACG
      thisacg = ccgR(:,excitedcells(j),excitedcells(j));
      axh = AxesInsetBars(a,.2,[.5 .5 .5],tR',thisacg);
      axhpos = get(axh,'Position');
      set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)],'XTickLabel',[]);
    end
    SquareSubplot(size(excitedcells,1)+1,size(excitedcells,1)+1);
    bar(tR,ccgR(:,refcell,refcell),'k');
    xlim([-0.04 0.04]);
    xlabel('Reference Cell ACG');
    set(fg,'Position',position);
    if ~isempty(db)
      DBConnect
      DBUse(db);
      try
	DBAddFigure(fg,['Hpc-Amygdala'],[xml ' shank ' int2str(refshankcell(1)) ' cell ' int2str(refshankcell(2)) ' CCG-20ms-Sorted Excitatory Connections'],'','','MonoSynConv.m');
      catch
      end
    end
    if ~strcmp(folder,'none')
      print(fg,[folder '/' xml '-shank' int2str(refshankcell(1)) '-cell' int2str(refshankcell(2)) '-CCG-20ms-Sorted Excitatory Connections'] ,'-dpng');
    end
    close(fg)
  end
%
  for i=1:size(FinalInhCellList,1)
    refcell=FinalInhCellList(i);
    refshankcell=completeIndex(completeIndex(:,3)==refcell,1:2);
    inhibitedcells = FinalInhMonoSynList(FinalInhMonoSynList(:,1)==refcell,2);
    fg=figure;
    suptitle([xml ' - Putative Inhibitory Cell : shank ' int2str(completeIndex(completeIndex(:,3)==refcell,1)) ' cell ' int2str(completeIndex(completeIndex(:,3)==refcell,2)) ' (ID : ' int2str(refcell) ')']);
    for j=1:size(inhibitedcells,1)
      a=SquareSubplot(size(inhibitedcells,1)+1,j);
      bar(tR, ccgR(:,refcell,inhibitedcells(j)));
      hold on;
      plot(tR,Pred(:,refcell,inhibitedcells(j)),'g');
      % Plot upper and lower boundaries
      PlotHVLines(Bounds(1,refcell,inhibitedcells(j)),'h','r--');
      PlotHVLines(Bounds(2,refcell,inhibitedcells(j)),'h','r--');
      % Plot signif increased bins in red
      exc=ccgR(:,refcell,inhibitedcells(j));
      exc(exc<Bounds(1,refcell,inhibitedcells(j))|~window)=0;
      bar(tR,exc,'r');
      % Plot signif lower bins in blue
      inh=ccgR(:,refcell,inhibitedcells(j));
      inh(inh>Bounds(2,refcell,inhibitedcells(j))|~window)=0;
      bar(tR,inh,'c');
      xlim([-0.02 0.02]);
      PlotColorIntervals([0.00075 0.00375],'LightGrey','v');
      xlabel(['shank ' int2str(completeIndex(completeIndex(:,3)==inhibitedcells(j),1)) ' cell ' int2str(completeIndex(completeIndex(:,3)==inhibitedcells(j),2)) ' (ID:' int2str(inhibitedcells(j)) ')' ]);
      % Plot an inset with the inhibitedcell ACG
      thisacg = ccgR(:,inhibitedcells(j),inhibitedcells(j));
      axh = AxesInsetBars(a,.2,[.5 .5 .5],tR',thisacg);
      axhpos = get(axh,'Position');
      set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)],'XTickLabel',[]);
    end
    SquareSubplot(size(inhibitedcells,1)+1,size(inhibitedcells,1)+1);
    bar(tR,ccgR(:,refcell,refcell),'k');
    xlim([-0.04 0.04]);
    xlabel('Reference Cell ACG');
    set(gcf,'Position',position);
    if ~isempty(db)
      DBConnect
      DBUse(db);
      try
	DBAddFigure(fg,['Hpc-Amygdala'],[xml ' shank ' int2str(refshankcell(1)) ' cell ' int2str(refshankcell(2)) ' CCG-20ms-Sorted Inhibitory Connections'],'','','MonoSynConv.m');
      catch
      end
    end
    if ~strcmp(folder,'none')
      print(fg,[folder '/' xml '-shank' int2str(refshankcell(1)) '-cell' int2str(refshankcell(2)) '-CCG-20ms-Sorted Inhibitory Connections'] ,'-dpng');
    end
    close (fg)
  end
end

FinalExcCellListID=FinalExcCellList;
FinalExcMonoSynID=FinalExcMonoSynList;
if ~isempty(FinalExcCellListID)
  FinalExcCellList=ConvertIDtoCell(FinalExcCellListID,completeIndex);
  FinalExcMonoSyn=ConvertIDtoCell(FinalExcMonoSynList,completeIndex);
else
  FinalExcCellList=[];
  FinalExcMonoSyn=[];
end
FinalInhCellListID=FinalInhCellList;
FinalInhMonoSynID=FinalInhMonoSynList;
if ~isempty(FinalInhCellListID)
  FinalInhCellList=ConvertIDtoCell(FinalInhCellListID,completeIndex);
  FinalInhMonoSyn=ConvertIDtoCell(FinalInhMonoSynList,completeIndex);
else
  FinalInhCellList=[];
  FinalInhMonoSyn=[];
end
Idx=completeIndex;
save([xml '-MonoSynConvClick.mat'],'Bounds','FinalExcCellListID','FinalExcCellList','FinalInhCellListID','FinalInhCellList','FinalExcMonoSynID','FinalExcMonoSyn','FinalInhMonoSynID','FinalInhMonoSyn','oops','Idx','ccgR','tR','Pred');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subplotclick(obj,ev) %when an axes is clicked

%  figobj = get(obj,'Parent');
axdata = get(obj,'UserData');

global SpuriousExcMonoSynList;
global SpuriousInhMonoSynList;

clr = get(obj,'Color');
if sum(clr == [1 1 1])==3%if white (ie synapse), set to pink (bad), remember as bad
    set(obj,'Color',[1 .75 .75])
    if strcmp(axdata.type,'exc')
      SpuriousExcMonoSynList(end+1,:)=axdata.cell;
    end
    if strcmp(axdata.type,'inh')
      SpuriousInhMonoSynList(end+1,:)=axdata.cell;
    end
elseif sum(clr == [1 .75 .75])==3%if pink, set to white, unremember as bad
    set(obj,'Color',[1 1 1])
    if strcmp(axdata.type,'exc')
      SpuriousExcMonoSynList(SpuriousExcMonoSynList(:,1)==axdata.cell(1)&SpuriousExcMonoSynList(:,2)==axdata.cell(2),:)=[];
    end
    if strcmp(axdata.type,'inh')
      SpuriousInhMonoSynList(SpuriousInhMonoSynList(:,1)==axdata.cell(1)&SpuriousInhMonoSynList(:,2)==axdata.cell(2),:)=[];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CloseFigureFcnExc(obj,ev)
%
global db
global xml
% cell=get(obj,'UserData')
figobj=get(obj,'Parent');
cell = guidata(obj);
%
if ~isempty(db)
  DBConnect
  DBUse(db);
  try
    DBAddFigure(gcf,['Hpc-Amygdala'],[xml ' shank ' int2str(cell(1)) ' cell ' int2str(cell(2)) ' CCG-20ms-Raw Excitatory Connections'],'','','MonoSynConvClick.m');
  catch
    disp ('Could not save figure to database')
  end
end
delete(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CloseFigureFcnInh(obj,ev)
%
global db
global xml
% cell=get(obj,'UserData')
cell = guidata(obj);
%
if ~isempty(db)
  DBConnect
  DBUse(db);
  try
    DBAddFigure(gcf,['Hpc-Amygdala'],[xml ' shank ' int2str(cell(1)) ' cell ' int2str(cell(2)) ' CCG-20ms-Raw Inhibitory Connections'],'','','MonoSynConvClick.m');
  catch
    disp ('Could not save figure to database')
  end
end
delete(obj);
