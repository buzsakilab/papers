function ExplainedVariance_PrePostCorr_Ultimate(struct,varargin)

% ExplainedVariance_PrePostCorr_Ultimate - Calculates and plots many variations of distrbutions of correlation changes between pre and post-sleep.
%   
%  USAGE
%
%    ExplainedVariance_PrepostCorr_Ultimate (struct,<options>)
%
%    struct		structure name (ex : 'BLA')
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%    EVtype		'cross' (cross-structure : default), 'intra' (intrastructure) 
%    celltype		'all' (default), 'pyr'
%    sleep		'sws' (default), 'rem','run'
%    =========================================================================
%
%  OUTPUT
%
%    Plots
%  
%  NOTE
%
%  SEE
%
%    See also : binspikes, corrcoef, corr
%
% January 2016, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
EVtype='cross';
celltype='all';
sleep='sws';

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
	switch(varargin{i}),
	  case 'EVtype',
		EVtype = lower(varargin{i+1});
	  case 'celltype'
		celltype = lower(varargin{i+1});
	  case 'sleep'
		sleep = lower(varargin{i+1});
	  otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

load('/media/Data-01/All-Rats/sessionindexing.mat')
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');

EVREV=[];
sess=[];
Allcorr.pre.isup=[];
Allcorr.pre.isdown=[];
Allcorr.pre.nomod=[];
Allcorr.pre.all=[];
Allcorr.run.isup=[];
Allcorr.run.isdown=[];
Allcorr.run.nomod=[];
Allcorr.run.all=[];
Allcorr.post.isup=[];
Allcorr.post.isdown=[];
Allcorr.post.nomod=[];
Allcorr.post.all=[];

AllpVal.pre.isup=[];
AllpVal.pre.isdown=[];
AllpVal.pre.nomod=[];
AllpVal.pre.all=[];
AllpVal.run.isup=[];
AllpVal.run.isdown=[];
AllpVal.run.nomod=[];
AllpVal.run.all=[];
AllpVal.post.isup=[];
AllpVal.post.isdown=[];
AllpVal.post.nomod=[];
AllpVal.post.all=[];


for i=1:length(xmlpath)
  xml=xmlpath{i}(end-14:end-1);
  cd(xmlpath{i});
  ratsess=ratsessionindex(strcmp(xmlpath{i},xmlpath),1:2); %[rat session]
  
  if strcmp(EVtype,'cross') & exist([xml '-ExplainedVariance-' struct '.mat'])==2 & exist([xml '.rip.evt'])==2
    load([xml '-ExplainedVariance-' struct '.mat'])
    treatsess=1;
    if isfield(CorrMatrix.pre,'run')
      treatrun=1;
    else
      treatrun=0;
    end
  elseif strcmp(EVtype,'intra') & exist([xml '-ExplainedVariance-' struct '-IntraStructure.mat'])==2 & exist([xml '.rip.evt'])==2
    load([xml '-ExplainedVariance-' struct '-IntraStructure.mat']);
    treatsess=1;
  else
    treatsess=0;
  end
   
  if treatsess
    if strcmp(sleep,'sws')
      CorrM.pre=CorrMatrix.pre.sws.all.Rall;
      CorrM.post=CorrMatrix.post.sws.all.Rall;
      CorrM.run=CorrMatrix.run.all.Rall;
      PvalM.pre=PvalMatrix.pre.sws.all.Rall;
      PvalM.post=PvalMatrix.post.sws.all.Rall;
      PvalM.run=PvalMatrix.run.all.Rall;
    elseif strcmp(sleep,'rem')
      CorrM.pre=CorrMatrix.pre.rem.all.Rall;
      CorrM.post=CorrMatrix.post.rem.all.Rall;
      CorrM.run=CorrMatrix.run.all.Rall;
      PvalM.pre=PvalMatrix.pre.rem.all.Rall;
      PvalM.post=PvalMatrix.post.rem.all.Rall;
      PvalM.run=PvalMatrix.run.all.Rall;
    elseif strcmp(sleep,'run')
      if treatrun
	CorrM.pre=CorrMatrix.pre.run.all.Rall;
	CorrM.post=CorrMatrix.post.run.all.Rall;
	CorrM.run=CorrMatrix.run.all.Rall;
	PvalM.pre=PvalMatrix.pre.run.all.Rall;
	PvalM.post=PvalMatrix.post.run.all.Rall;
	PvalM.run=PvalMatrix.run.all.Rall;
      else
	CorrM.pre=[];
	CorrM.post=[];
	CorrM.run=[];
	PvalM.pre=[];
	PvalM.post=[];
	PvalM.run=[];
      end
    end
  
    if strcmp(EVtype,'cross')
      if strcmp(celltype,'pyr');
	STispyr=STtype==1;
	STselectcells.up=STispyr&is.up;
	STselectcells.down=STispyr&is.down;
	STselectcells.nomod=STispyr&~is.mod;
	STselectcells.all=STispyr;
	HPispyr=HPtype==1;
	HPselectcells=HPispyr;
      elseif strcmp(celltype,'all')
	STselectcells.up=is.up;
	STselectcells.down=is.down;
	STselectcells.nomod=~is.mod;
	STselectcells.all=':';
	HPselectcells=':';
      end
    elseif strcmp(EVtype,'intra');
      if strcmp(celltype,'pyr')
	STselectcells.up=STispyr&is.up;
	STselectcells.down=STispyr&is.down;
	STselectcells.nomod=STispyr&~is.mod;
	STselectcells.all=STispyr;
      elseif strcmp(celltype,'all')
	STselectcells.up=is.up;
	STselectcells.down=is.down;
	STselectcells.nomod=~is.mod;
	STselectcells.all=':';
      end
    end
      
    if strcmp(EVtype,'cross')
      if ~isempty(CorrM.pre)
	CorrM.pre_isup=CorrM.pre(HPselectcells,STselectcells.up);
	Allcorr.pre.isup=[Allcorr.pre.isup;CorrM.pre_isup(:)];
	CorrM.pre_isdown=CorrM.pre(HPselectcells,STselectcells.down);
	Allcorr.pre.isdown=[Allcorr.pre.isdown;CorrM.pre_isdown(:)];
	CorrM.pre_nomod=CorrM.pre(HPselectcells,STselectcells.nomod);
	Allcorr.pre.nomod=[Allcorr.pre.nomod; CorrM.pre_nomod(:)];
	CorrM.pre_all=CorrM.pre(HPselectcells,STselectcells.all);
	Allcorr.pre.all=[Allcorr.pre.all; CorrM.pre_all(:)];
	
	CorrM.run_isup=CorrM.run(HPselectcells,STselectcells.up);
	Allcorr.run.isup=[Allcorr.run.isup;CorrM.run_isup(:)];
	CorrM.run_isdown=CorrM.run(HPselectcells,STselectcells.down);
	Allcorr.run.isdown=[Allcorr.run.isdown;CorrM.run_isdown(:)];
	CorrM.run_nomod=CorrM.run(HPselectcells,STselectcells.nomod);
	Allcorr.run.nomod=[Allcorr.run.nomod; CorrM.run_nomod(:)];
	CorrM.run_all=CorrM.run(HPselectcells,STselectcells.all);
	Allcorr.run.all=[Allcorr.run.all;CorrM.run_all(:)];
	
	CorrM.post_isup=CorrM.post(HPselectcells,STselectcells.up);
	Allcorr.post.isup=[Allcorr.post.isup;CorrM.post_isup(:)];
	CorrM.post_isdown=CorrM.post(HPselectcells,STselectcells.down);
	Allcorr.post.isdown=[Allcorr.post.isdown;CorrM.post_isdown(:)];
	CorrM.post_nomod=CorrM.post(HPselectcells,STselectcells.nomod);
	Allcorr.post.nomod=[Allcorr.post.nomod; CorrM.post_nomod(:)];
	CorrM.post_all=CorrM.post(HPselectcells,STselectcells.all);
	Allcorr.post.all=[Allcorr.post.all;CorrM.post_all(:)];
	
	%%
	PvalM.pre_isup=PvalM.pre(HPselectcells,STselectcells.up);
	AllpVal.pre.isup=[AllpVal.pre.isup;PvalM.pre_isup(:)];
	PvalM.pre_isdown=PvalM.pre(HPselectcells,STselectcells.down);
	AllpVal.pre.isdown=[AllpVal.pre.isdown;PvalM.pre_isdown(:)];
	PvalM.pre_nomod=PvalM.pre(HPselectcells,STselectcells.nomod);
	AllpVal.pre.nomod=[AllpVal.pre.nomod; PvalM.pre_nomod(:)];
	PvalM.pre_all=PvalM.pre(HPselectcells,STselectcells.all);
	AllpVal.pre.all=[AllpVal.pre.all; PvalM.pre_all(:)];
	
	PvalM.run_isup=PvalM.run(HPselectcells,STselectcells.up);
	AllpVal.run.isup=[AllpVal.run.isup;PvalM.run_isup(:)];
	PvalM.run_isdown=PvalM.run(HPselectcells,STselectcells.down);
	AllpVal.run.isdown=[AllpVal.run.isdown;PvalM.run_isdown(:)];
	PvalM.run_nomod=PvalM.run(HPselectcells,STselectcells.nomod);
	AllpVal.run.nomod=[AllpVal.run.nomod; PvalM.run_nomod(:)];
	PvalM.run_all=PvalM.run(HPselectcells,STselectcells.all);
	AllpVal.run.all=[AllpVal.run.all;PvalM.run_all(:)];
	
	PvalM.post_isup=PvalM.post(HPselectcells,STselectcells.up);
	AllpVal.post.isup=[AllpVal.post.isup;PvalM.post_isup(:)];
	PvalM.post_isdown=PvalM.post(HPselectcells,STselectcells.down);
	AllpVal.post.isdown=[AllpVal.post.isdown;PvalM.post_isdown(:)];
	PvalM.post_nomod=PvalM.post(HPselectcells,STselectcells.nomod);
	AllpVal.post.nomod=[AllpVal.post.nomod; PvalM.post_nomod(:)];
	PvalM.post_all=PvalM.post(HPselectcells,STselectcells.all);
	AllpVal.post.all=[AllpVal.post.all;PvalM.post_all(:)];
      end
    elseif strcmp(EVtype,'intra')     
      if sum(STselectcells.up)>2
	CorrM.pre_isup=CorrM.pre(STselectcells.up,STselectcells.up);
	Allcorr.pre.isup=[Allcorr.pre.isup;CorrM.pre_isup(logical(triu((ones(size(CorrM.pre_isup))),1)))];
	CorrM.pre_isdown=CorrM.pre(STselectcells.down,STselectcells.down);
	Allcorr.pre.isdown=[Allcorr.pre.isdown;CorrM.pre_isdown(logical(triu((ones(size(CorrM.pre_isdown))),1)))];
	CorrM.pre_nomod=CorrM.pre(STselectcells.nomod,STselectcells.nomod);
	Allcorr.pre.nomod=[Allcorr.pre.nomod; CorrM.pre_nomod(logical(triu((ones(size(CorrM.pre_nomod))),1)))];
	CorrM.pre_all=CorrM.pre(STselectcells.all,STselectcells.all);
	Allcorr.pre.all=[Allcorr.pre.all; CorrM.pre_all(logical(triu((ones(size(CorrM.pre_all))),1)))];
	
	CorrM.run_isup=CorrM.run(STselectcells.up,STselectcells.up);
	Allcorr.run.isup=[Allcorr.run.isup;CorrM.run_isup(logical(triu((ones(size(CorrM.run_isup))),1)))];
	CorrM.run_isdown=CorrM.run(STselectcells.down,STselectcells.down);
	Allcorr.run.isdown=[Allcorr.run.isdown;CorrM.run_isdown(logical(triu((ones(size( CorrM.run_isdown))),1)))];
	CorrM.run_nomod=CorrM.run(STselectcells.nomod,STselectcells.nomod);
	Allcorr.run.nomod=[Allcorr.run.nomod; CorrM.run_nomod(logical(triu((ones(size(CorrM.run_nomod))),1)))];
	CorrM.run_all=CorrM.run(STselectcells.all,STselectcells.all);
	Allcorr.run.all=[Allcorr.run.all;CorrM.run_all(logical(triu((ones(size(CorrM.run_all))),1)))];
	
	CorrM.post_isup=CorrM.post(STselectcells.up,STselectcells.up);
	Allcorr.post.isup=[Allcorr.post.isup;CorrM.post_isup(logical(triu((ones(size(CorrM.post_isup))),1)))];
	CorrM.post_isdown=CorrM.post(STselectcells.down,STselectcells.down);
	Allcorr.post.isdown=[Allcorr.post.isdown;CorrM.post_isdown(logical(triu((ones(size(CorrM.post_isdown))),1)))];
	CorrM.post_nomod=CorrM.post(STselectcells.nomod,STselectcells.nomod);
	Allcorr.post.nomod=[Allcorr.post.nomod; CorrM.post_nomod(logical(triu((ones(size(CorrM.post_nomod))),1)))];
	CorrM.post_all=CorrM.post(STselectcells.all,STselectcells.all);
	Allcorr.post.all=[Allcorr.post.all;CorrM.post_all(logical(triu((ones(size(CorrM.post_all))),1)))];
	
	%%
	PvalM.pre_isup=PvalM.pre(STselectcells.up,STselectcells.up);
	AllpVal.pre.isup=[AllpVal.pre.isup;PvalM.pre_isup(logical(triu((ones(size(PvalM.pre_isup))),1)))];
	PvalM.pre_isdown=PvalM.pre(STselectcells.down,STselectcells.down);
	AllpVal.pre.isdown=[AllpVal.pre.isdown;PvalM.pre_isdown(logical(triu((ones(size(PvalM.pre_isdown))),1)))];
	PvalM.pre_nomod=PvalM.pre(STselectcells.nomod,STselectcells.nomod);
	AllpVal.pre.nomod=[AllpVal.pre.nomod; PvalM.pre_nomod(logical(triu((ones(size(PvalM.pre_nomod))),1)))];
	PvalM.pre_all=PvalM.pre(STselectcells.all,STselectcells.all);
	AllpVal.pre.all=[AllpVal.pre.all; PvalM.pre_all(logical(triu((ones(size(PvalM.pre_all))),1)))];
	
	PvalM.run_isup=PvalM.run(STselectcells.up,STselectcells.up);
	AllpVal.run.isup=[AllpVal.run.isup;PvalM.run_isup(logical(triu((ones(size(PvalM.run_isup))),1)))];
	PvalM.run_isdown=PvalM.run(STselectcells.down,STselectcells.down);
	AllpVal.run.isdown=[AllpVal.run.isdown;PvalM.run_isdown(logical(triu((ones(size(PvalM.run_isdown))),1)))];
	PvalM.run_nomod=PvalM.run(STselectcells.nomod,STselectcells.nomod);
	AllpVal.run.nomod=[AllpVal.run.nomod; PvalM.run_nomod(logical(triu((ones(size( PvalM.run_nomod))),1)))];
	PvalM.run_all=PvalM.run(STselectcells.all,STselectcells.all);
	AllpVal.run.all=[AllpVal.run.all;PvalM.run_all(logical(triu((ones(size(PvalM.run_all))),1)))];
	
	PvalM.post_isup=PvalM.post(STselectcells.up,STselectcells.up);
	AllpVal.post.isup=[AllpVal.post.isup;PvalM.post_isup(logical(triu((ones(size(PvalM.post_isup))),1)))];
	PvalM.post_isdown=PvalM.post(STselectcells.down,STselectcells.down);
	AllpVal.post.isdown=[AllpVal.post.isdown;PvalM.post_isdown(logical(triu((ones(size(PvalM.post_isdown))),1)))];
	PvalM.post_nomod=PvalM.post(STselectcells.nomod,STselectcells.nomod);
	AllpVal.post.nomod=[AllpVal.post.nomod; PvalM.post_nomod(logical(triu((ones(size(PvalM.post_nomod))),1)))];
	PvalM.post_all=PvalM.post(STselectcells.all,STselectcells.all);
	AllpVal.post.all=[AllpVal.post.all;PvalM.post_all(logical(triu((ones(size(PvalM.post_all))),1)))];	
      end
    end
  end
end

prepostDiff.isup=Allcorr.post.isup-Allcorr.pre.isup;
prepostDiff.isdown=Allcorr.post.isdown-Allcorr.pre.isdown;
prepostDiff.nomod=Allcorr.post.nomod-Allcorr.pre.nomod;
prepostDiff.all=Allcorr.post.all-Allcorr.pre.all;


% Test for normality
%  [h,p]=kstest(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01))
%  [h,p]=kstest(prepostDiff.isup(AllpVal.run.isup>0.01))
%  [h,p]=kstest(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01))
%  [h,p]=kstest(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01))
%  [h,p]=kstest(prepostDiff.nomod(AllpVal.run.nomod>0.01))
%  [h,p]=kstest(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01))
%  [h,p]=kstest(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01))
%  [h,p]=kstest(prepostDiff.isdown(AllpVal.run.isdown>0.01))
%  [h,p]=kstest(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01))

%%ANOVA and shit.
ripplegr=[ones(length(prepostDiff.isup),1);ones(length(prepostDiff.nomod),1)*2;ones(length(prepostDiff.isdown),1)*3];
anovadiff=[prepostDiff.isup;prepostDiff.nomod;prepostDiff.isdown];
runpval=[AllpVal.run.isup;AllpVal.run.nomod;AllpVal.run.isdown];
runcorr=[Allcorr.run.isup;Allcorr.run.nomod;Allcorr.run.isdown];
% 1=upmod, 2=nomod, 3=downmod
issigplus=runpval<0.01&runcorr>0;
issigneg=runpval<0.01&runcorr<0;
isnosig=runpval>=0.01;

percentsigplus=sum(issigplus)/(sum(issigplus)+sum(isnosig)+sum(issigneg))*100
percentsigneg=sum(issigneg)/(sum(issigplus)+sum(isnosig)+sum(issigneg))*100
percentnosig=sum(isnosig)/(sum(issigplus)+sum(isnosig)+sum(issigneg))*100

anovagr=zeros(length(runpval),1);

anovagr(ripplegr==1&issigplus)=1;
mean(anovadiff(ripplegr==1))
anovagr(ripplegr==1&isnosig)=2;
anovagr(ripplegr==1&issigneg)=3;
anovagr(ripplegr==2&issigplus)=4;
anovagr(ripplegr==2&isnosig)=5;
anovagr(ripplegr==2&issigneg)=6;
anovagr(ripplegr==3&issigplus)=7;
anovagr(ripplegr==3&isnosig)=8;
anovagr(ripplegr==3&issigneg)=9;

anovadiff(anovagr==0)=[];
anovagr(anovagr==0)=[];


[p,anovatab,stats]=anova1(anovadiff,anovagr,'Alpha',0.001);
%  [p,anovatab,stats]=kruskalwallis(anovadiff,anovagr);
[comp,means]=multcompare(stats,'Alpha',0.01);

%Equality of variance test
p=vartestn(anovadiff,anovagr);


%%%%%%%%%%%%
if strcmp(celltype,'all')
  celltypetitle='Pyr+Int';
elseif strcmp(celltype,'pyr')
  celltypetitle='PyrOnly';
end

if strcmp(EVtype,'cross')
  evtypetitle='Cross-Structure';
elseif strcmp(EVtype,'intra')
  evtypetitle='Intra-Structure';
end

anovagr=[];

f1=figure('Position',[110 67 1841 936]); %PRE-POST correlation distributions - Per ripple-mod type and RUN correlation significance.

subplot(4,4,1);
[h,bins]=hist(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01),100);
bar(bins,h,'FaceColor',rgb('Crimson'),'EdgeColor','none');
xlabel(['upmod cells, sig RUN corr + n=' num2str(sum(h))]);
PlotHVLines(0,'r');
mean_isup_runcorrpos=nanmean(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01));
sem_isup_runcorrpos=nansem(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01));

subplot(4,4,4);hold on;
CumulativePlot(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.1]);PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');
subplot(4,4,13);hold on;
CumulativePlot(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.1]);PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');

subplot(4,4,5);
h2=hist(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01),bins);
bar(bins,h2,'FaceColor',rgb('FireBrick'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['upmod cells, sig RUN corr - n=' num2str(sum(h2))]);
mean_isup_runcorrneg=nanmean(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01));
sem_isup_runcorrneg=nansem(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01));

subplot(4,4,8);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01),'normalize','on','color','FireBrick','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');
subplot(4,4,13);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01),'normalize','on','color','FireBrick','newfig','off');

subplot(4,4,9);
h3=hist(prepostDiff.isup(AllpVal.run.isup>0.01),bins);
bar(bins,h3,'FaceColor',rgb('DarkGrey'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['upmod cells, NON sig RUN corr n=' num2str(sum(h3))]);
mean_isup_runcorrno=nanmean(prepostDiff.isup(AllpVal.run.isup>0.01));
sem_isup_runcorrno=nansem(prepostDiff.isup(AllpVal.run.isup>0.01));

subplot(4,4,12);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isup(AllpVal.run.isup>0.01),'normalize','on','color','DarkGrey','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');
subplot(4,4,13);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isup(AllpVal.run.isup>0.01),'normalize','on','color','DarkGrey','newfig','off');

subplot(4,4,2);
h3=hist(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01),bins);
bar(bins,h3,'FaceColor',rgb('SteelBlue'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['downmod cells, sig RUN corr + n=' num2str(sum(h3))]);
mean_isdown_runcorrpos=nanmean(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01));
sem_isdown_runcorrpos=nansem(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01));

subplot(4,4,4);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01),'normalize','on','color','SteelBlue','newfig','off');
subplot(4,4,14);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01),'normalize','on','color','SteelBlue','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');

subplot(4,4,6);
h4=hist(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01),bins);
bar(bins,h4,'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['downmod cells, sig RUN corr - n=' num2str(sum(h4))]);
mean_isdown_runcorrneg=nanmean(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01));
sem_isdown_runcorrneg=nansem(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01));

subplot(4,4,8);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01),'normalize','on','color','DodgerBlue','newfig','off');
subplot(4,4,14);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01),'normalize','on','color','DodgerBlue','newfig','off');

subplot(4,4,10);
h5=hist(prepostDiff.isdown(AllpVal.run.isdown>0.01),bins);
bar(bins,h5,'FaceColor',rgb('LightBlue'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['downmod cells, Non sig RUN corr n=' num2str(sum(h5))]);
mean_isdown_runcorrno=nanmean(prepostDiff.isdown(AllpVal.run.isdown>0.01));
sem_isdown_runcorrno=nansem(prepostDiff.isdown(AllpVal.run.isdown>0.01));

subplot(4,4,12);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(AllpVal.run.isdown>0.01),'normalize','on','color','LightBlue','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');
subplot(4,4,14);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.isdown(AllpVal.run.isdown>0.01),'normalize','on','color','LightBlue','newfig','off');

subplot(4,4,3);
h5=hist(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01),bins);
bar(bins,h5,'FaceColor',rgb('SeaGreen'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['nomod cells, sig RUN corr + n=' int2str(sum(h5))]);
mean_nomod_runcorrpos=nanmean(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01));
sem_nomod_runcorrpos=nansem(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01));

subplot(4,4,4);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01),'normalize','on','color','SeaGreen','newfig','off');
subplot(4,4,15);hold on;;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01),'normalize','on','color','SeaGreen','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');

subplot(4,4,7);
h6=hist(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01),bins);
bar(bins,h6,'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['nomod cells, sig RUN corr - n=' int2str(sum(h6))]);
mean_nomod_runcorrneg=nanmean(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01));
sem_nomod_runcorrneg=nansem(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01));

subplot(4,4,8);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01),'normalize','on','color','MediumSeaGreen','newfig','off');
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');
subplot(4,4,15);hold on;;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01),'normalize','on','color','MediumSeaGreen','newfig','off');
subplot(4,4,11);
h7=hist(prepostDiff.nomod(AllpVal.run.nomod>0.01),bins);
bar(bins,h7,'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');
PlotHVLines(0,'r');
xlabel(['nomod cells, Non sig RUN corr n=' int2str(sum(h7))]);
mean_nomod_runcorrno=nanmean(prepostDiff.nomod(AllpVal.run.nomod>0.01));
sem_nomod_runcorrno=nansem(prepostDiff.nomod(AllpVal.run.nomod>0.01));

subplot(4,4,12);hold on;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(AllpVal.run.nomod>0.01),'normalize','on','color','DarkSeaGreen','newfig','off');
subplot(4,4,15);hold on;;xlim([-0.1 0.1]);
CumulativePlot(prepostDiff.nomod(AllpVal.run.nomod>0.01),'normalize','on','color','DarkSeaGreen','newfig','off');
suptitle(['Prepost correlation difference distributions - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);


%%%%

f1b=figure;hold on;
barwitherr([sem_isup_runcorrpos 0 0 0 0 0 0 0 0 0 0],[mean_isup_runcorrpos 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');hold on;
barwitherr([0 sem_isup_runcorrno 0 0 0 0 0 0 0 0 0],[0 mean_isup_runcorrno 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkGrey'),'EdgeColor','none');hold on;
barwitherr([0 0 sem_isup_runcorrneg 0 0 0 0 0 0 0 0],[0 0 mean_isup_runcorrneg 0 0 0 0 0 0 0 0],'FaceColor',rgb('FireBrick'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 sem_nomod_runcorrpos 0 0 0 0 0 0],[0 0 0 0 mean_nomod_runcorrpos 0 0 0 0 0 0],'FaceColor',rgb('SeaGreen'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 0 sem_nomod_runcorrno 0 0 0 0 0],[0 0 0 0 0 mean_nomod_runcorrno 0 0 0 0 0],'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 0 0 sem_nomod_runcorrneg 0 0 0 0],[0 0 0 0 0 0 mean_nomod_runcorrneg 0 0 0 0],'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 0 0 0 0 sem_isdown_runcorrpos 0 0],[0 0 0 0 0 0 0 0 mean_isdown_runcorrpos 0 0],'FaceColor',rgb('SteelBlue'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 0 0 0 0 0 sem_isdown_runcorrno 0],[0 0 0 0 0 0 0 0 0 mean_isdown_runcorrno 0],'FaceColor',rgb('LightBlue'),'EdgeColor','none');hold on;
barwitherr([0 0 0 0 0 0 0 0 0 0 sem_isdown_runcorrneg],[0 0 0 0 0 0 0 0 0 0 mean_isdown_runcorrneg],'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');
xlabel(['anova p=' num2str(p)]);





% Prepare fucking boxplot
f1c=figure;
data=[prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01);
prepostDiff.isup(AllpVal.run.isup>0.01);
prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01);
prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01);
prepostDiff.nomod(AllpVal.run.nomod>0.01);
prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01);    
prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01);
prepostDiff.isdown(AllpVal.run.isdown>0.01);
prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01)];
groups=[ones(length(prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01)),1);
ones(length(prepostDiff.isup(AllpVal.run.isup>0.01)),1)*2;
ones(length(prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01)),1)*3;
ones(length(prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01)),1)*4;
ones(length(prepostDiff.nomod(AllpVal.run.nomod>0.01)),1)*5;
ones(length(prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01)),1)*6;    
ones(length(prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01)),1)*7;
ones(length(prepostDiff.isdown(AllpVal.run.isdown>0.01)),1)*8;
ones(length(prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01)),1)*9];
boxplot(data,groups); %It's awful.

f1c=figure;
violin({prepostDiff.isup(Allcorr.run.isup>0&AllpVal.run.isup<0.01) prepostDiff.isup(AllpVal.run.isup>0.01) prepostDiff.isup(Allcorr.run.isup<0&AllpVal.run.isup<0.01) prepostDiff.nomod(Allcorr.run.nomod>0&AllpVal.run.nomod<0.01) prepostDiff.nomod(AllpVal.run.nomod>0.01) prepostDiff.nomod(Allcorr.run.nomod<0&AllpVal.run.nomod<0.01) prepostDiff.isdown(Allcorr.run.isdown>0&AllpVal.run.isdown<0.01) prepostDiff.isdown(AllpVal.run.isdown>0.01) prepostDiff.isdown(Allcorr.run.isdown<0&AllpVal.run.isdown<0.01)}) 

means=[mean_isup_runcorrpos;mean_isup_runcorrno;mean_isup_runcorrneg;mean_nomod_runcorrpos;mean_nomod_runcorrno;mean_nomod_runcorrneg;mean_isdown_runcorrpos;mean_isdown_runcorrno;mean_isdown_runcorrneg]
%%%%%%%%%%%%%%%
f2=figure('Position',[321 182 1441 775]);	%Pre/Post Correlation difference for different run correlations.
suptitle(['Prepost corr diff distr - RUN sig corr -'  struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);
subplot(2,3,1);
[h,bins]=hist(prepostDiff.all(Allcorr.run.all>0&AllpVal.run.all<0.01),100);
bar(bins,h,'FaceColor',rgb('DarkOrange'),'EdgeColor','none');
xlabel('Sig. Positive RUN corr');
PlotHVLines(0,'r');
subplot(2,3,3);
[h,bins]=hist(prepostDiff.all(Allcorr.run.all<0&AllpVal.run.all<0.01),bins);
bar(bins,h,'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
xlabel('Sig. Negative RUN corr');
PlotHVLines(0,'r');
subplot(2,3,2);
[h,bins]=hist(prepostDiff.all(AllpVal.run.all>0.01),bins);
bar(bins,h,'FaceColor',rgb('DarkGrey'),'EdgeColor','none');
xlabel('NON Sig.RUN corr');
PlotHVLines(0,'r');
subplot(2,3,4:6)
CumulativePlot(prepostDiff.all(Allcorr.run.all>0&AllpVal.run.all<0.01),'normalize','on','color','DarkOrange','newfig','off');
hold on;
CumulativePlot(prepostDiff.all(Allcorr.run.all<0&AllpVal.run.all<0.01),'normalize','on','color','MidnightBlue','newfig','off');
CumulativePlot(prepostDiff.all(AllpVal.run.all>0.01),'normalize','on','color','DarkGrey','newfig','off');
legend('Sig. Pos. Run corr','Sig. Neg. run corr','Non sig run corr','Location','NorthWest');
PlotHVLines(0,'v','k:');
PlotHVLines(0.5,'h','k:');
xlabel('Post-Pre Sleep Correlation Difference');


%%%%%%%%%%%%%%%
f3=figure('Position',[70 78 1625 877]);
suptitle(['Prepost corr. difference distribs - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);
subplot(3,4,1);
[h,bins]=hist(prepostDiff.isup,100);
bar(bins,h,'FaceColor',rgb('IndianRed'),'EdgeColor','none');
xlabel('All upmod')
PlotHVLines(0,'r');
subplot(3,4,2);
h2=hist(prepostDiff.isdown,bins);
bar(bins,h2,'FaceColor',rgb('RoyalBlue'),'EdgeColor','none');
xlabel('All downmod');
PlotHVLines(0,'r');
subplot(3,4,3);
h3=hist(prepostDiff.nomod,bins);
bar(bins,h3,'FaceColor',rgb('DarkGray'),'EdgeColor','none');
xlabel('All nomod');
PlotHVLines(0,'r');
subplot(3,4,4);hold on;
CumulativePlot(prepostDiff.isup,'normalize','on','color','IndianRed','newfig','off');
CumulativePlot(prepostDiff.isdown,'normalize','on','color','RoyalBlue','newfig','off');
CumulativePlot(prepostDiff.nomod,'normalize','on','color','DarkGray','newfig','off');
PlotHVLines(0,'v','k:');
PlotHVLines(0.5,'h','k:');
xlim([-0.1 0.1]);

subplot(3,4,5);
[h,bins]=hist(prepostDiff.isup(AllpVal.run.isup<0.01),100);
bar(bins,h,'FaceColor',rgb('Crimson'),'EdgeColor','none');
xlabel('upmod, sig. corr during RUN')
PlotHVLines(0,'r');
subplot(3,4,6);
h2=hist(prepostDiff.isdown(AllpVal.run.isdown<0.01),bins);
bar(bins,h2,'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
xlabel('downmod, sig corr during run');
PlotHVLines(0,'r');
subplot(3,4,7);
h3=hist(prepostDiff.nomod(AllpVal.run.nomod<0.01),bins);
bar(bins,h3,'FaceColor',rgb('DimGray'),'EdgeColor','none');
xlabel('nomod, sig corr run');
PlotHVLines(0,'r');
subplot(3,4,8);hold on;
CumulativePlot(prepostDiff.isup(AllpVal.run.isup<0.01),'normalize','on','color','Crimson','newfig','off');
CumulativePlot(prepostDiff.isdown(AllpVal.run.isdown<0.01),'normalize','on','color','MidnightBlue','newfig','off');
CumulativePlot(prepostDiff.nomod(AllpVal.run.nomod<0.01),'normalize','on','color','DimGray','newfig','off');
xlim([-0.1 0.1]);
PlotHVLines(0,'v','k:');
PlotHVLines(0.5,'h','k:');

subplot(3,4,9);hold on;
CumulativePlot(prepostDiff.isup,'normalize','on','color','IndianRed','newfig','off');
CumulativePlot(prepostDiff.isup(AllpVal.run.isup<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.1]);
PlotHVLines(0,'v','k:');
PlotHVLines(0.5,'h','k:');
subplot(3,4,10);hold on;
CumulativePlot(prepostDiff.isdown,'normalize','on','color','RoyalBlue','newfig','off');
CumulativePlot(prepostDiff.isdown(AllpVal.run.isdown<0.01),'normalize','on','color','MidnightBlue','newfig','off');
xlim([-0.1 0.1]);
PlotHVLines(0,'v','k:');
PlotHVLines(0.5,'h','k:');
subplot(3,4,11);hold on;
CumulativePlot(prepostDiff.nomod,'normalize','on','color','DarkGray','newfig','off');
CumulativePlot(prepostDiff.nomod(AllpVal.run.nomod<0.01),'normalize','on','color','DimGray','newfig','off');
xlim([-0.1 0.1]);
PlotHVLines(0,'v','k:');PlotHVLines(0.5,'h','k:');


%%%%%%%%%%%%%%%
f4=figure('Position',[1 88 1786 876]);
suptitle(['Pre/Post Correlations correlation + distributions - All corr - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);

subplot(3,6,7)%%
hold on;
plot(Allcorr.pre.isup,Allcorr.post.isup,'.');
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.isup,Allcorr.post.isup,'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h(1),'color','r')
plot([-0.05 0.1],[-0.05 0.1],'k:');
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('UpMod cells, all corr');
subplot(3,6,1)%%1
[h,bins]=hist(Allcorr.pre.isup,100);
bar(bins,h,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
subplot(3,6,8)%%
h2=hist(Allcorr.post.isup,bins);
barh(bins,h2,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,2)%%
bar(bins,h2-h,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
xlabel('post-pre hist');
subplot(3,6,13:14);hold on;
CumulativePlot(Allcorr.pre.isup,'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.isup,'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

subplot(3,6,9)
plot(Allcorr.pre.isdown,Allcorr.post.isdown,'.');
hold on
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.isdown,Allcorr.post.isdown,'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h,'color','r');
plot([-0.05 0.1],[-0.05 0.1],'k:')
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('DownMod cells, all corr');
subplot(3,6,3)
h3=hist(Allcorr.pre.isdown,bins);
bar(bins,h3,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
ylim([0 500]);
subplot(3,6,10)
h4=hist(Allcorr.post.isdown,bins);
barh(bins,h4,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,4)
bar(bins,h4-h3,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
xlim([-0.05 0.1]);
hold on; PlotHVLines(0,'r');
subplot(3,6,15:16);hold on;
CumulativePlot(Allcorr.pre.isdown,'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.isdown,'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

subplot(3,6,11)
plot(Allcorr.pre.nomod,Allcorr.post.nomod,'.');
hold on
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.nomod,Allcorr.post.nomod,'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h,'color','r');
plot([-0.05 0.1],[-0.05 0.1],'k:');
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('NoMod cells, all corr');
subplot(3,6,5)
h5=hist(Allcorr.pre.nomod,bins);
bar(bins,h5,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
xlim([-0.05 0.1]);
ylim([0 7000]);
hold on; PlotHVLines(0,'r');
subplot(3,6,12)
h6=hist(Allcorr.post.nomod,bins);
barh(bins,h6,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
xlim([0 7000]);
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,6);
bar(bins,h6-h5,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
xlim([-0.05 0.1]);
hold on; PlotHVLines(0,'r');
subplot(3,6,17:18);hold on;
CumulativePlot(Allcorr.pre.nomod,'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.nomod,'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5=figure('Position',[1 88 1786 876]);
suptitle(['Pre/Post Correlations correlation + distributions - Sig. run corr only - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);

subplot(3,6,7)%%
hold on;
plot(Allcorr.pre.isup(AllpVal.run.isup<0.01),Allcorr.post.isup(AllpVal.run.isup<0.01),'.');
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.isup(AllpVal.run.isup<0.01),Allcorr.post.isup(AllpVal.run.isup<0.01),'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h(1),'color','r')
plot([-0.05 0.1],[-0.05 0.1],'k:');
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('UpMod cells, all corr');
subplot(3,6,1)%%1
[h,bins]=hist(Allcorr.pre.isup(AllpVal.run.isup<0.01),100);
bar(bins,h,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
subplot(3,6,8)%%
h2=hist(Allcorr.post.isup(AllpVal.run.isup<0.01),bins);
barh(bins,h2,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,2)%%
bar(bins,h2-h,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
xlabel('post-pre hist');
subplot(3,6,13:14);hold on;
CumulativePlot(Allcorr.pre.isup(AllpVal.run.isup<0.01),'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.isup(AllpVal.run.isup<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

subplot(3,6,9)
plot(Allcorr.pre.isdown(AllpVal.run.isdown<0.01),Allcorr.post.isdown(AllpVal.run.isdown<0.01),'.');
hold on
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.isdown(AllpVal.run.isdown<0.01),Allcorr.post.isdown(AllpVal.run.isdown<0.01),'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h,'color','r');
plot([-0.05 0.1],[-0.05 0.1],'k:')
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('DownMod cells, all corr');
subplot(3,6,3)
h3=hist(Allcorr.pre.isdown(AllpVal.run.isdown<0.01),bins);
bar(bins,h3,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
hold on; PlotHVLines(0,'r');
xlim([-0.05 0.1]);
subplot(3,6,10)
h4=hist(Allcorr.post.isdown(AllpVal.run.isdown<0.01),bins);
barh(bins,h4,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,4)
bar(bins,h4-h3,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
xlim([-0.05 0.1]);
hold on; PlotHVLines(0,'r');
subplot(3,6,15:16);hold on;
CumulativePlot(Allcorr.pre.isdown(AllpVal.run.isdown<0.01),'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.isdown(AllpVal.run.isdown<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

subplot(3,6,11)
plot(Allcorr.pre.nomod(AllpVal.run.nomod<0.01),Allcorr.post.nomod(AllpVal.run.nomod<0.01),'.');
hold on
xlim([-0.05 0.1]);
ylim([-0.05 0.1]);
[r1,p1]=corrcoef(Allcorr.pre.nomod(AllpVal.run.nomod<0.01),Allcorr.post.nomod(AllpVal.run.nomod<0.01),'rows','complete');
r1=r1(1,2);p1=p1(1,2);
h=lsline;
set(h,'color','r');
plot([-0.05 0.1],[-0.05 0.1],'k:');
title(['r:' num2str(r1) ' p:' num2str(p1)]);
xlabel('NoMod cells, all corr');
subplot(3,6,5)
h5=hist(Allcorr.pre.nomod(AllpVal.run.nomod<0.01),bins);
bar(bins,h5,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('DarkGreen'));
xlim([-0.05 0.1]);
hold on; PlotHVLines(0,'r');
subplot(3,6,12)
h6=hist(Allcorr.post.nomod(AllpVal.run.nomod<0.01),bins);
barh(bins,h6,'EdgeColor','none','BarWidth',1.1,'FaceColor',rgb('Crimson'));
hold on; PlotHVLines(0,'h','r');
ylim([-0.05 0.1]);
subplot(3,6,6);
bar(bins,h6-h5,'EdgeColor','none','BarWidth',1.1,'FaceColor','k');
xlim([-0.05 0.1]);
hold on; PlotHVLines(0,'r');
subplot(3,6,17:18);hold on;
CumulativePlot(Allcorr.pre.nomod(AllpVal.run.nomod<0.01),'normalize','on','color','DarkGreen','newfig','off');
CumulativePlot(Allcorr.post.nomod(AllpVal.run.nomod<0.01),'normalize','on','color','Crimson','newfig','off');
xlim([-0.1 0.15]);
PlotHVLines(0,'v',':');
PlotHVLines(0.5,'h',':');

%%%%%%%%%%%%%
f6=figure('Position',[225 415 1360 358]);  % Correlation between RUN correlation and pre/post increase/decrease in corr for ripple modulated cells
suptitle(['RUN corr vs prepost changes - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);
subplot(1,3,1)
plot(Allcorr.run.isup,prepostDiff.isup,'.');
xlim([-0.2 0.5]);ylim([-0.12 0.12]);
h=lsline;
set(h,'color','r');
xlabel('RUN corr');
ylabel('Pre/post Corr difference');
title('Up-mod cells');
subplot(1,3,2)
plot(Allcorr.run.isdown,prepostDiff.isdown,'.');
h=lsline;
set(h,'color','r');
xlim([-0.2 0.5]);ylim([-0.12 0.12]);
xlabel('RUN corr');
ylabel('Prepost Corr diff');
title('Down-mod cells');
subplot(1,3,3)
plot(Allcorr.run.nomod,prepostDiff.nomod,'.');
h=lsline;
set(h,'color','r');
xlim([-0.2 0.5]);ylim([-0.12 0.12]);
xlabel('RUN corr');
ylabel('Prepost Corr diff');
title('No-mod cells');

f7=figure('Position',[172 129 1360 850]);
suptitle(['PRE corr vs POST corr - ' struct ' - ' evtypetitle ' - ' celltypetitle ' - ' sleep]);
subplot(2,3,1);hold on;
plot(Allcorr.pre.all,Allcorr.post.all,'.','color',rgb('Black'),'MarkerSize',15);
h=lsline;
set(h,'color','r');
plot([-0.15 0.3],[-0.15 0.3],'k:');
xlim([-0.15 0.3]);ylim([-0.15 0.3]);
xlabel('PRE corr');
ylabel('POST corr');
title('All corr');
subplot(2,3,2);hold on;
plot(Allcorr.pre.all(AllpVal.run.all<0.01),Allcorr.post.all(AllpVal.run.all<0.01),'.','color',rgb('Crimson'),'MarkerSize',15);
h=lsline;
set(h,'color','r');
plot([-0.15 0.3],[-0.15 0.3],'k:');
xlim([-0.15 0.3]);ylim([-0.15 0.3]);
xlabel('PRE corr');
ylabel('POST corr');
title('Sig. run corr');
subplot(2,3,3);hold on;
plot(Allcorr.pre.all(AllpVal.run.all>0.01),Allcorr.post.all(AllpVal.run.all>0.01),'.','color',rgb('Gray'),'MarkerSize',15);
h=lsline;
set(h,'color','r');
plot([-0.15 0.3],[-0.15 0.3],'k:');
xlim([-0.15 0.3]);ylim([-0.15 0.3]);
xlabel('PRE corr');
ylabel('POST corr');
title('Non sig run corr');
subplot(2,3,5);hold on;
plot(Allcorr.pre.all(Allcorr.run.all>0 & AllpVal.run.all<0.01),Allcorr.post.all(Allcorr.run.all>0 & AllpVal.run.all<0.01),'.','color',rgb('DarkOrange'),'MarkerSize',15);
h=lsline;
set(h,'color','r');
plot([-0.15 0.3],[-0.15 0.3],'k:');
xlim([-0.15 0.3]);ylim([-0.15 0.3]);
xlabel('PRE corr');
ylabel('POST corr');
title('Sig + run corr');
subplot(2,3,6);hold on;
plot(Allcorr.pre.all(Allcorr.run.all<0 & AllpVal.run.all<0.01),Allcorr.post.all(Allcorr.run.all<0 & AllpVal.run.all<0.01),'.','color',rgb('MidnightBlue'),'MarkerSize',15);
h=lsline;
set(h,'color','r');
plot([-0.15 0.3],[-0.15 0.3],'k:');
xlim([-0.15 0.3]);ylim([-0.15 0.3]);
xlabel('PRE corr');
ylabel('POST corr');
title('Sig - run corr');


savefig=input('Save Figures?','s')
if strcmp(savefig,'yes')
  cd('/media/Data-01/All-Rats/AllRats-RunPrePost-PairwiseCorrelations');
  figure(f1)
  plot2svg(['Figure1-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure1-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure1-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f1b)
  plot2svg(['Figure1b-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure1b-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure1b-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f2)
  plot2svg(['Figure2-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure2-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure2-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f3)
  plot2svg(['Figure3-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure3-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure3-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f4)
  plot2svg(['Figure4-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure4-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure4-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f5)
  plot2svg(['Figure5-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure5-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure5-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f6)
  plot2svg(['Figure6-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure6-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure6-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
  figure(f7)
  plot2svg(['Figure7-' struct '-' celltypetitle '-' evtypetitle '-' sleep '.svg'],gcf);
  saveas(gcf,['Figure7-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'png');
  saveas(gcf,['Figure7-' struct '-' celltypetitle '-' evtypetitle '-' sleep],'pdf');
end