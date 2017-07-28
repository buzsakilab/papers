function [Allcorr,Allpval,Index,EV,REV,PairContrib,Minusindex,BLA,HPC,CorrDiff] = ExplainedVariance_PooledCorr(struct,threshold,varargin)

% Explained Variance_PooledCorr - Calculates global EV and REV, pair and cell contributions. Plots RUN vs sleep correlation changes for various contribution groups.
%  
%  USAGE
%
%    [AllCorr,Index,RippleMod] = ExplainedVariance_PooledCorr (struct,<options>)
%
%    struct		structure name (ex : 'BLA')
%    threshold		threshold for contributing and anticontributing cell PAIRS
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%    EVtype		'cross' (cross-structure : default), 'intra' (intrastructure)
%    thresh2		threshold for mean contribution (per cell) Default = 0,000003
%    =========================================================================
%
%  OUTPUT
    
%      Allcorr              matrices of correlation values (triplets) pre/RUN/post for sws, rem sleep, run and ripples only.
%      Allpval              associated pvalues for Allcorr
%      EV                   global EV
%      REV                  gloval REV
%      Paircontrib          individual pair contributions
%      MinusIndex           Index for pair contributions CorrMatrix
%      BLA                  different groups, lists of BLA cells
%      HPC                  different groups, lists of HPC cells
%      CorrDiff             difference in correlations between pre and post rem/sws/run/ripples
%  
%  NOTE
%     
%      This program calculates variable 'AllRats-EVcontribs.mat'
%  
%  SEE
%
% March 2016, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
EVtype='cross';
thresh2=0.000003;
%  celltype='all';
%  sleep='sws';

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
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

load('/media/Data-01/All-Rats/sessionindexing.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/PoissonRippleMod.mat');
load('/media/Data-01/All-Rats/AllRats-StatesFR.mat');

Allcorr.sws=[];
Allcorr.rem=[];
Allcorr.run=[];
Allcorr.rip=[];
Allpval.sws=[];
Allpval.rem=[];
Allpval.run=[];
Allpval.rip=[];

Index.Cells=[];
Index.Type=[];
Index.RippleMod=[];
Index.FR=[];

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
    Allcorr.sws=[Allcorr.sws;CorrMatrix.pre.sws.all.Rall(:) CorrMatrix.run.all.Rall(:) CorrMatrix.post.sws.all.Rall(:)];
    Allcorr.rem=[Allcorr.rem;CorrMatrix.pre.rem.all.Rall(:) CorrMatrix.run.all.Rall(:) CorrMatrix.post.rem.all.Rall(:)];
    Allcorr.rip=[Allcorr.rip;CorrMatrix.pre.rip_in.all.Rall(:) CorrMatrix.run.all.Rall(:) CorrMatrix.post.rip_in.all.Rall(:)];
    if treatrun
      Allcorr.run=[Allcorr.run;CorrMatrix.pre.run.all.Rall(:) CorrMatrix.run.all.Rall(:) CorrMatrix.post.run.all.Rall(:)];
      Allpval.run=[Allpval.run;PvalMatrix.pre.run.all.Rall(:) PvalMatrix.run.all.Rall(:) PvalMatrix.post.run.all.Rall(:)];
    else
      Allcorr.run=[Allcorr.run;NaN(length(CorrMatrix.pre.sws.all.Rall(:)),3)];
      Allpval.run=[Allpval.run;NaN(length(PvalMatrix.pre.sws.all.Rall(:)),3)];
    end
    Allpval.sws=[Allpval.sws;PvalMatrix.pre.sws.all.Rall(:) PvalMatrix.run.all.Rall(:) PvalMatrix.post.sws.all.Rall(:)];
    Allpval.rem=[Allpval.rem;PvalMatrix.pre.rem.all.Rall(:) PvalMatrix.run.all.Rall(:) PvalMatrix.post.rem.all.Rall(:)];
    Allpval.rip=[Allpval.rip;PvalMatrix.pre.rip_in.all.Rall(:) PvalMatrix.run.all.Rall(:) PvalMatrix.post.rip_in.all.Rall(:)];
    STindex=STindex(:,1:2);
    HPindex=HPindex(:,1:2);
    Index.Cells=[Index.Cells;[repmat(ratsess,length(CorrMatrix.pre.sws.all.Rall(:)),1) repmat(HPindex(:,1:2),size(STindex,1),1) reshape(repmat(STindex(:)',size(HPindex,1),1),[],2)]];
    Index.Type=[Index.Type;[repmat(HPtype,size(STindex,1),1) reshape(repmat(STtype',size(HPindex,1),1),[],1)]];
    
    STripplemod=zeros(size(STindex,1),1);
    STripplemod(is.up)=1;
    STripplemod(is.down)=2;    
    Index.RippleMod=[Index.RippleMod; reshape(repmat(STripplemod',size(HPindex,1),1),[],1)];
    
    STfiringrate=StatesFR(ismember(StatesFR(:,1:4),[repmat(ratsess,size(STindex,1),1) STindex],'rows'),10);
    HPfiringrate=StatesFR(ismember(StatesFR(:,1:4),[repmat(ratsess,size(HPindex,1),1) HPindex],'rows'),10);
    Index.FR=[Index.FR;[repmat(HPfiringrate,size(STindex,1),1) reshape(repmat(STfiringrate',size(HPindex,1),1),[],1)]];    
  end
end

bothpyr=sum(Index.Type,2)==2;
nonan.rem=sum(isnan(Allcorr.rem),2)==0;
nonan.sws=sum(isnan(Allcorr.sws),2)==0;
nonan.run=sum(isnan(Allcorr.run),2)==0;

r_runpost.rem.all=corrcoef(Allcorr.rem(nonan.rem,3),Allcorr.rem(nonan.rem,2));r_runpost.rem.all=r_runpost.rem.all(1,2);
r_runpre.rem.all=corrcoef(Allcorr.rem(nonan.rem,1),Allcorr.rem(nonan.rem,2));r_runpre.rem.all=r_runpre.rem.all(1,2);
r_prepost.rem.all=corrcoef(Allcorr.rem(nonan.rem,1),Allcorr.rem(nonan.rem,3));r_prepost.rem.all=r_prepost.rem.all(1,2);

r_runpost.run.all=corrcoef(Allcorr.run(nonan.run,3),Allcorr.run(nonan.run,2));r_runpost.run.all=r_runpost.run.all(1,2);
r_runpre.run.all=corrcoef(Allcorr.run(nonan.run,1),Allcorr.run(nonan.run,2));r_runpre.run.all=r_runpre.run.all(1,2);
r_prepost.run.all=corrcoef(Allcorr.run(nonan.run,1),Allcorr.run(nonan.run,3));r_prepost.run.all=r_prepost.run.all(1,2);

r_runpost.sws.all=corrcoef(Allcorr.sws(nonan.sws,3),Allcorr.sws(nonan.sws,2));r_runpost.sws.all=r_runpost.sws.all(1,2);
r_runpre.sws.all=corrcoef(Allcorr.sws(nonan.sws,1),Allcorr.sws(nonan.sws,2));r_runpre.sws.all=r_runpre.sws.all(1,2);
r_prepost.sws.all=corrcoef(Allcorr.sws(nonan.sws,1),Allcorr.sws(nonan.sws,3));r_prepost.sws.all=r_prepost.sws.all(1,2);

r_runpost.rem.pyr=corrcoef(Allcorr.rem(nonan.rem&bothpyr,3),Allcorr.rem(nonan.rem&bothpyr,2));r_runpost.rem.pyr=r_runpost.rem.pyr(1,2);
r_runpre.rem.pyr=corrcoef(Allcorr.rem(nonan.rem&bothpyr,1),Allcorr.rem(nonan.rem&bothpyr,2));r_runpre.rem.pyr=r_runpre.rem.pyr(1,2);
r_prepost.rem.pyr=corrcoef(Allcorr.rem(nonan.rem&bothpyr,1),Allcorr.rem(nonan.rem&bothpyr,3));r_prepost.rem.pyr=r_prepost.rem.pyr(1,2);

r_runpost.sws.pyr=corrcoef(Allcorr.sws(nonan.sws&bothpyr,3),Allcorr.sws(nonan.sws&bothpyr,2));r_runpost.sws.pyr=r_runpost.sws.pyr(1,2);
r_runpre.sws.pyr=corrcoef(Allcorr.sws(nonan.sws&bothpyr,1),Allcorr.sws(nonan.sws&bothpyr,2));r_runpre.sws.pyr=r_runpre.sws.pyr(1,2);
r_prepost.sws.pyr=corrcoef(Allcorr.sws(nonan.sws&bothpyr,1),Allcorr.sws(nonan.sws&bothpyr,3));r_prepost.sws.pyr=r_prepost.sws.pyr(1,2);

r_runpost.run.pyr=corrcoef(Allcorr.run(nonan.run&bothpyr,3),Allcorr.run(nonan.run&bothpyr,2));r_runpost.run.pyr=r_runpost.run.pyr(1,2);
r_runpre.run.pyr=corrcoef(Allcorr.run(nonan.run&bothpyr,1),Allcorr.run(nonan.run&bothpyr,2));r_runpre.run.pyr=r_runpre.run.pyr(1,2);
r_prepost.run.pyr=corrcoef(Allcorr.run(nonan.run&bothpyr,1),Allcorr.run(nonan.run&bothpyr,3));r_prepost.run.pyr=r_prepost.run.pyr(1,2);


EV.rem.pyr=((r_runpost.rem.pyr-r_runpre.rem.pyr*r_prepost.rem.pyr)/sqrt((1-r_runpre.rem.pyr^2)*(1-r_prepost.rem.pyr^2)))^2;
REV.rem.pyr=((r_runpre.rem.pyr-r_runpost.rem.pyr*r_prepost.rem.pyr)/sqrt((1-r_runpost.rem.pyr^2)*(1-r_prepost.rem.pyr^2)))^2;
EV.rem.all=((r_runpost.rem.all-r_runpre.rem.all*r_prepost.rem.all)/sqrt((1-r_runpre.rem.all^2)*(1-r_prepost.rem.all^2)))^2;
REV.rem.all=((r_runpre.rem.all-r_runpost.rem.all*r_prepost.rem.all)/sqrt((1-r_runpost.rem.all^2)*(1-r_prepost.rem.all^2)))^2;

EV.sws.pyr=((r_runpost.sws.pyr-r_runpre.sws.pyr*r_prepost.sws.pyr)/sqrt((1-r_runpre.sws.pyr^2)*(1-r_prepost.sws.pyr^2)))^2;
REV.sws.pyr=((r_runpre.sws.pyr-r_runpost.sws.pyr*r_prepost.sws.pyr)/sqrt((1-r_runpost.sws.pyr^2)*(1-r_prepost.sws.pyr^2)))^2;
EV.sws.all=((r_runpost.sws.all-r_runpre.sws.all*r_prepost.sws.all)/sqrt((1-r_runpre.sws.all^2)*(1-r_prepost.sws.all^2)))^2;
REV.sws.all=((r_runpre.sws.all-r_runpost.sws.all*r_prepost.sws.all)/sqrt((1-r_runpost.sws.all^2)*(1-r_prepost.sws.all^2)))^2;

EV.run.pyr=((r_runpost.run.pyr-r_runpre.run.pyr*r_prepost.run.pyr)/sqrt((1-r_runpre.run.pyr^2)*(1-r_prepost.run.pyr^2)))^2;
REV.run.pyr=((r_runpre.run.pyr-r_runpost.run.pyr*r_prepost.run.pyr)/sqrt((1-r_runpost.run.pyr^2)*(1-r_prepost.run.pyr^2)))^2;
EV.run.all=((r_runpost.run.all-r_runpre.run.all*r_prepost.run.all)/sqrt((1-r_runpre.run.all^2)*(1-r_prepost.run.all^2)))^2;
REV.run.all=((r_runpre.run.all-r_runpost.run.all*r_prepost.run.all)/sqrt((1-r_runpost.run.all^2)*(1-r_prepost.run.all^2)))^2;


figure;
subplot(2,1,1)
bar([EV.sws.all REV.sws.all 0 EV.rem.all REV.rem.all 0 EV.run.all REV.run.all]);
subplot(2,1,2)
xlabel('AllCells');
bar([EV.sws.pyr REV.sws.pyr 0 EV.rem.pyr REV.rem.pyr 0 EV.run.pyr REV.run.pyr]);
xlabel('Pyr only')
suptitle('Global EV/REV : SWS - REM - RUN');

EV.minus=[];
REV.minus=[];
Minusindex=[];

% Contributions calc for pyr-pyr pairs only
for i=1:length(Allcorr.sws)
  keep=ones(length(Allcorr.sws),1);
  if nonan.sws(i) & bothpyr(i)
    keep(i)=0; % changing the cell to remove to 0 in keep (= keep all but ith cell)
 
    r_runpost.minus=corrcoef(Allcorr.sws(nonan.sws&bothpyr&keep,3),Allcorr.sws(nonan.sws&bothpyr&keep,2));r_runpost.minus=r_runpost.minus(1,2);
    r_runpre.minus=corrcoef(Allcorr.sws(nonan.sws&bothpyr&keep,1),Allcorr.sws(nonan.sws&bothpyr&keep,2));r_runpre.minus=r_runpre.minus(1,2);
    r_prepost.minus=corrcoef(Allcorr.sws(nonan.sws&bothpyr&keep,1),Allcorr.sws(nonan.sws&bothpyr&keep,3));r_prepost.minus=r_prepost.minus(1,2);

    EV.minus=[EV.minus ; ((r_runpost.minus-r_runpre.minus*r_prepost.minus)/sqrt((1-r_runpre.minus^2)*(1-r_prepost.minus^2)))^2];
    REV.minus=[REV.minus ; ((r_runpre.minus-r_runpost.minus*r_prepost.minus)/sqrt((1-r_runpost.minus^2)*(1-r_prepost.minus^2)))^2];
    
    Minusindex=[Minusindex;Index.Cells(i,:)];
  end
end

EV.diffs=EV.minus-EV.sws.pyr;

Index.EVdiffs=NaN(length(Index.Cells),1);
Index.EVdiffs(ismember(Index.Cells,Minusindex,'rows'))=EV.diffs;

%%% Threshold pair contrib
PairContrib.Tplus=Minusindex(EV.diffs<-threshold,:);
PairContrib.Tminus=Minusindex(EV.diffs>threshold,:);
PairContrib.Taverageplus=Minusindex(EV.diffs<=0&EV.diffs>-threshold,:);
PairContrib.Taverageminus=Minusindex(EV.diffs>0&EV.diffs<threshold,:);
PairContrib.TBLAplus=unique(PairContrib.Tplus(:,[1 2 5 6]),'rows');
PairContrib.TBLAminus=unique(PairContrib.Tminus(:,[1 2 5 6]),'rows');
PairContrib.THPCplus=unique(PairContrib.Tplus(:,[1:4]),'rows');
PairContrib.THPCminus=unique(PairContrib.Tminus(:,[1:4]),'rows');

%%%% Quartiles pair contrib
contribqu=quantile(EV.diffs,[0.25 0.5 0.75]);
PairContrib.Qplus=Minusindex(EV.diffs<contribqu(1),:);
PairContrib.Qminus=Minusindex(EV.diffs>contribqu(3),:);
PairContrib.Qaverageplus=Minusindex(EV.diffs>contribqu(1)&EV.diffs<contribqu(2),:);
PairContrib.Qaverageminus=Minusindex(EV.diffs<contribqu(3)&EV.diffs>contribqu(2),:);
PairContrib.QBLAplus=unique(PairContrib.Qplus(:,[1 2 5 6]),'rows');
PairContrib.QBLAminus=unique(PairContrib.Qminus(:,[1 2 5 6]),'rows');
PairContrib.QHPCplus=unique(PairContrib.Qplus(:,[1:4]),'rows');
PairContrib.QHPCminus=unique(PairContrib.Qminus(:,[1:4]),'rows');

%%%% x percent paircontrib !!! variable named 5p but various percentiles hard-coded
contrib5p=prctile(EV.diffs,[2.5 97.5]) %percentile hard-coded here
PairContrib.p5plus=Minusindex(EV.diffs<contrib5p(1),:);
PairContrib.p5minus=Minusindex(EV.diffs>contrib5p(2),:);
PairContrib.p5average=Minusindex(EV.diffs>contrib5p(1)&EV.diffs<contrib5p(2),:);

figure;
[h,bins]=hist(EV.diffs,5000);
bar(bins,h,'FaceColor',rgb('Gray'),'EdgeColor','none');
hold on;
PlotHVLines(contribqu,'v','r');
PlotHVLines([threshold*(-1);threshold],'v','k');
PlotHVLines(contrib5p,'v','g')
xlabel('Change in EV');
ylabel('# cell pairs');

%PerCellContrib
BLAindex=unique(Minusindex(:,[1 2 5 6]),'rows');
HPCindex=unique(Minusindex(:,1:4),'rows');

BLA.all=BLAindex;
HPC.all=HPCindex;

for i=1:length(BLAindex)
  BLAmeancontribs(i)=mean(EV.diffs(ismember(Minusindex(:,[1 2 5 6]),BLAindex(i,:),'rows')));
end

for i=1:length(HPCindex)
  HPCmeancontribs(i)=mean(EV.diffs(ismember(Minusindex(:,1:4),HPCindex(i,:),'rows')));
end

%%% Firing rate GRoups
BLAFR=StatesFR(ismember(StatesFR(:,1:4),BLAindex,'rows'),10);
HPCFR=StatesFR(ismember(StatesFR(:,1:4),HPCindex,'rows'),10);
% BLA
blaqu=quantile(BLAFR,[0.25 0.5 0.75]);
blahighFRQ=BLAFR>blaqu(3);
blalowFRQ=BLAFR<blaqu(1);
blamedhighFRQ=BLAFR>blaqu(1)&BLAFR<blaqu(2);
blamedlowFRQ=BLAFR>blaqu(2)&BLAFR<blaqu(3);
BLA.lowFR=BLAindex(blalowFRQ,:);
BLA.highFR=BLAindex(blahighFRQ,:);
BLA.medlowFR=BLAindex(blamedlowFRQ,:);
BLA.medhighFR=BLAindex(blamedhighFRQ,:);
% Hpc
hpcqu=quantile(HPCFR,[0.25 0.5 0.75]);
hpchighFRQ=HPCFR>hpcqu(3);
hpclowFRQ=HPCFR<hpcqu(1);
hpcmedhighFRQ=HPCFR>hpcqu(1)&HPCFR<hpcqu(2);
hpcmedlowFRQ=HPCFR>hpcqu(2)&HPCFR<hpcqu(3);
HPC.lowFR=HPCindex(hpclowFRQ,:);
HPC.highFR=HPCindex(hpchighFRQ,:);
HPC.medlowFR=HPCindex(hpcmedlowFRQ,:);
HPC.medhighFR=HPCindex(hpcmedhighFRQ,:);

%%% Contribution groups-BLA
BLA.q=quantile(BLAmeancontribs,[0.25 0.5 0.75]);
BLAhighQuant=BLAmeancontribs<BLA.q(1);
BLAlowQuant=BLAmeancontribs>BLA.q(3);
BLAMedhighQuant=BLAmeancontribs>BLA.q(1)&BLAmeancontribs<BLA.q(2);
BLAMedlowQuant=BLAmeancontribs>BLA.q(2)&BLAmeancontribs<BLA.q(3);
% quantiles
BLA.highQ=BLAindex(BLAhighQuant,:);
BLA.lowQ=BLAindex(BLAlowQuant,:);
BLA.medhighQ=BLAindex(BLAMedhighQuant,:);
BLA.medlowQ=BLAindex(BLAMedlowQuant,:);
% threshold
BLA.highT=BLAindex(BLAmeancontribs<-thresh2,:);
BLA.lowT=BLAindex(BLAmeancontribs>thresh2,:);
BLA.averageT=BLAindex(BLAmeancontribs>-thresh2&BLAmeancontribs<thresh2,:);

%%% Contribution groups-HPC
HPC.q=quantile(HPCmeancontribs,[0.25 0.5 0.75]);
HPChighQuant=HPCmeancontribs<HPC.q(1);
HPClowQuant=HPCmeancontribs>HPC.q(3);
HPCMedhighQuant=HPCmeancontribs>HPC.q(1)&HPCmeancontribs<HPC.q(2);
HPCMedlowQuant=HPCmeancontribs>HPC.q(2)&HPCmeancontribs<HPC.q(3);
% quantiles
HPC.highQ=HPCindex(HPChighQuant,:);
HPC.lowQ=HPCindex(HPClowQuant,:);
HPC.medhighQ=HPCindex(HPCMedhighQuant,:);
HPC.medlowQ=HPCindex(HPCMedlowQuant,:);
%threshold
HPC.highT=HPCindex(HPCmeancontribs<-thresh2,:);
HPC.lowT=HPCindex(HPCmeancontribs>thresh2,:);
HPC.averageT=HPCindex(HPCmeancontribs>-thresh2&HPCmeancontribs<thresh2,:);

figure;
subplot(1,2,1)
hist(BLAmeancontribs,200);
hold on;
PlotHVLines(BLA.q,'v','g:');
PlotHVLines([-thresh2 0 thresh2],'v','r:');
subplot(1,2,2)
hist(HPCmeancontribs,200);
PlotHVLines(HPC.q,'v','g');
PlotHVLines([-thresh2 0 thresh2],'v','r:');


cd('/media/Data-01/All-Rats/AllRats-EVcontributingCells');

CorrDiff.run=Allcorr.run(:,3)-Allcorr.run(:,1);
CorrDiff.rem=Allcorr.rem(:,3)-Allcorr.rem(:,1);
CorrDiff.sws=Allcorr.sws(:,3)-Allcorr.sws(:,1);
CorrDiff.rip=Allcorr.rip(:,3)-Allcorr.rip(:,1);

isNaN.run=isnan(CorrDiff.run);

%%% RUN/SWS correlation changes 
%%%%%%%%%%%%%%%%%%% threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure('Position',[302 543 1345 356]);
subplot(1,4,1)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Tplus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Tplus,'rows')),'.','Color',rgb('Crimson'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
lsline
%  subplot(1,4,2)
%  plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Taverageplus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Taverageplus,'rows')),'.','Color',rgb('Coral'),'MarkerSize',12);
%  xlim([-0.1 0.1]);ylim([-0.4 0.3]);
%  xlabel('Pre-Post diff. in corr. SWS');
%  ylabel('Pre-Post diff. in corr. testRUNS');
%  lsline
subplot(1,4,3)
plot(CorrDiff.sws(ismember(Index.Cells,[PairContrib.Taverageplus;PairContrib.Taverageminus],'rows')),CorrDiff.run(ismember(Index.Cells,[PairContrib.Taverageplus;PairContrib.Taverageminus],'rows')),'.','Color',rgb('DarkGrey'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
lsline
subplot(1,4,4)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Tminus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Tminus,'rows')),'.','Color',rgb('MidnightBlue'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
lsline
suptitle(['Pairwise Contributions to global EV - Threshold = ' num2str(threshold)]);

%%%%%%%%%%%%%%%%% x percentiles %%%%%%%%%%%%%%%%%%%%%%%%%%%
[r5plus,p5plus]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5plus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.p5plus,'rows')&~isNaN.run));
r5plus=r5plus(1,2)
p5plus=p5plus(1,2)
np5plus=length(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5plus,'rows')&~isNaN.run))
[r5av,p5av]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5average,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.p5average,'rows')&~isNaN.run));
r5av=r5av(1,2)
p5av=p5av(1,2)
np5av=length(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5average,'rows')&~isNaN.run))
[r5minus,p5minus]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5minus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.p5minus,'rows')&~isNaN.run));
r5minus=r5minus(1,2)
p5minus=p5minus(1,2)
np5minus=length(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5minus,'rows')&~isNaN.run))


figure('Position',[302 543 1345 356]);
subplot(1,4,1)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5plus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.p5plus,'rows')),'.','Color',rgb('Crimson'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(r5plus) ' p=' num2str(p5plus)]);
lsline
subplot(1,4,3)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5average,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.p5average,'rows')),'.','Color',rgb('DarkGrey'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(r5av) ' p=' num2str(p5av)]);
lsline
subplot(1,4,4)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.p5minus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.p5minus,'rows')),'.','Color',rgb('MidnightBlue'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(r5minus) ' p=' num2str(p5minus)]);
lsline
suptitle(['Pairwise Contributions to global EV 2.5 percent extremes : <' num2str(contrib5p(1)) ', >' num2str(contrib5p(2))]);

%%% RUN/SWS correlation changes
[rQplus,pQplus]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qplus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.Qplus,'rows')&~isNaN.run));
rQplus=rQplus(1,2)
pQplus=pQplus(1,2)
[rQavpl,pQavpl]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qaverageplus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.Qaverageplus,'rows')&~isNaN.run));
rQavpl=rQavpl(1,2)
pQavpl=pQavpl(1,2)
[rQavm,pQavm]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qaverageminus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.Qaverageminus,'rows')&~isNaN.run));
rQavm=rQavm(1,2)
pQavm=pQavm(1,2)
[rQminus,pQminus]=corrcoef(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qminus,'rows')&~isNaN.run),CorrDiff.run(ismember(Index.Cells,PairContrib.Qminus,'rows')&~isNaN.run));
rQminus=rQminus(1,2)
pQminus=pQminus(1,2)

figure('Position',[302 543 1345 356]);
subplot(1,4,1)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qplus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Qplus,'rows')),'.','Color',rgb('Crimson'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(rQplus) ' p=' num2str(pQplus)]);
lsline
subplot(1,4,2)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qaverageplus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Qaverageplus,'rows')),'.','Color',rgb('Coral'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(rQavpl) ' p=' num2str(pQavpl)]);
lsline
subplot(1,4,3)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qaverageminus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Qaverageminus,'rows')),'.','Color',rgb('RoyalBlue'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(rQavm) ' p=' num2str(pQavm)]);
lsline
subplot(1,4,4)
plot(CorrDiff.sws(ismember(Index.Cells,PairContrib.Qminus,'rows')),CorrDiff.run(ismember(Index.Cells,PairContrib.Qminus,'rows')),'.','Color',rgb('MidnightBlue'),'MarkerSize',12);
xlim([-0.1 0.1]);ylim([-0.4 0.3]);
xlabel('Pre-Post diff. in corr. SWS');
ylabel('Pre-Post diff. in corr. testRUNS');
title(['r=' num2str(rQminus) ' p=' num2str(pQminus)]);
lsline
suptitle('Pairwise Contributions to global EV - Quartiles');


