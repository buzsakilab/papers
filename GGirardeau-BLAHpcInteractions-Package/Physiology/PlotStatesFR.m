function [skewn,psignrank,med,Mean,Sem] = PlotStatesFR (structure,varargin)

%PlotStatesFR - Plots firing rate graphs for different states / structures
%
%  USAGE
%
%    PlotStatesFR (structure,<options>)
%
%    structure          Structure to plot from ex: 'BLA'
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'newfig'		Plot in a new figure (default = 'on')
%    =========================================================================
%
%  OUTPUT
%
%    Figures
%
%  NOTE
%
%
%  SEE
%
%    See also : StatesFR.m
%
% Gabrielle Girardeau, September 2014
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
newfig= 'on';

remfr=6;
swsfr=7;
wakefr=8;
position=[306 93 1138 881];

%  strucname=input('Structure?','s');

load('/media/Data-01/All-Rats/SpikeParameters.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/AllRats-StatesFR.mat');
load('/media/Data-01/All-Rats/Structures/structures.mat');

strucname=input('Structure?','s');

structure=eval(structure)

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
	  case 'newfig'
		  newfig = varargin{i+1};
	  otherwise,
		  error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

%Selection vectors
ccgType=SpikeParameters(:,end);
isccgPyr=ccgType==1;
isccgInt=ccgType==2;
isccgUndef=ccgType==0;
isStruc=ismember(StatesFR(:,1:3),structure,'rows');
finalType=finalType(:,end);
isfinalPyr=finalType==1;
isfinalInt=finalType==2;
isfinalUndef=finalType==0;

if strcmp(strucname,'CeCM')
  isccgPyr=zeros(7390,1);
  isccgInt=zeros(7390,1);
  isccgUndef=ones(7390,1);
end

%%% REM/wake
figure('Position',[487    68   928   913]);
subplot(2,2,3)
hold on
plot(log10(StatesFR(isStruc&isccgUndef,wakefr)),log10(StatesFR(isStruc&isccgUndef,remfr)),'.','color',rgb('DarkGray'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgPyr,wakefr)),log10(StatesFR(isStruc&isccgPyr,remfr)),'.','color',rgb('OrangeRed'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgInt,wakefr)),log10(StatesFR(isStruc&isccgInt,remfr)),'.','color',rgb('MediumBlue'),'MarkerSize',15);
set(gca,'XTick',[-3:2]);
set(gca,'YTick',[-3:2]);
xlim([-3 2])
ylim([-3 2])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
plot([-3 2],[-3 2],'r')
xlabel(['Wake']);
ylabel(['REM']);

subplot(2,2,1)
hold on
[h,bins]=hist(log10(StatesFR(isStruc&isccgUndef,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('DarkGray'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgPyr,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('OrangeRed'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgInt,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('MediumBlue'),'LineStyle','none');
set(gca,'XTick',[-3:2],'XTickLabel',[0.001 0.01 0.1 1 10 100]);
xlim([-3 2]);

subplot(2,2,4)
hold on
[h,bins]=hist(log10(StatesFR(isStruc&isccgUndef,remfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('DarkGray'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgPyr,remfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('OrangeRed'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgInt,remfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('MediumBlue'),'LineStyle','none');
set(gca,'YTick',[-3:2],'YTickLabel',[0.001 0.01 0.1 1 10 100]);
ylim([-3 2]);

suptitle(['All rats - State Firing Rates - REMWake' strucname]);
plot2svg(['Rats-All - ' strucname '-hist-REMWake.svg'],gcf);

%%% SWS/Wake
figure('Position',[487    68   928   913]);
subplot(2,2,3)
hold on
plot(log10(StatesFR(isStruc&isccgUndef,wakefr)),log10(StatesFR(isStruc&isccgUndef,swsfr)),'.','color',rgb('DarkGray'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgPyr,wakefr)),log10(StatesFR(isStruc&isccgPyr,swsfr)),'.','color',rgb('OrangeRed'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgInt,wakefr)),log10(StatesFR(isStruc&isccgInt,swsfr)),'.','color',rgb('MediumBlue'),'MarkerSize',15);
set(gca,'XTick',[-3:2]);
set(gca,'YTick',[-3:2]);
xlim([-3 2])
ylim([-3 2])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
plot([-3 2],[-3 2],'r')
xlabel(['Wake']);
ylabel(['REM']);

subplot(2,2,1)
hold on
[h,bins]=hist(log10(StatesFR(isStruc&isccgUndef,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('DarkGray'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgPyr,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('OrangeRed'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgInt,wakefr)),[-3:0.1:2]);
bar(bins,h,'FaceColor',rgb('MediumBlue'),'LineStyle','none');
set(gca,'XTick',[-3:2],'XTickLabel',[0.001 0.01 0.1 1 10 100]);
xlim([-3 2]);

subplot(2,2,4)
hold on
[h,bins]=hist(log10(StatesFR(isStruc&isccgUndef,swsfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('DarkGray'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgPyr,swsfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('OrangeRed'),'LineStyle','none');
[h,bins]=hist(log10(StatesFR(isStruc&isccgInt,swsfr)),[-3:0.1:2]);
barh(bins,h,'FaceColor',rgb('MediumBlue'),'LineStyle','none');
set(gca,'YTick',[-3:2],'YTickLabel',[0.001 0.01 0.1 1 10 100]);
ylim([-3 2]);

suptitle(['All rats - State Firing Rates - SWSwake' strucname]);
plot2svg(['Rats-All - ' strucname '-hist-SWSwake.svg'],gcf);


f0=figure('Position',position);

subplot(3,3,1)
hold on
plot(log10(StatesFR(isStruc&isccgUndef,wakefr)),log10(StatesFR(isStruc&isccgUndef,remfr)),'.','color',rgb('DarkGray'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgPyr,wakefr)),log10(StatesFR(isStruc&isccgPyr,remfr)),'.','color',rgb('OrangeRed'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgInt,wakefr)),log10(StatesFR(isStruc&isccgInt,remfr)),'.','color',rgb('MediumBlue'),'MarkerSize',15);
set(gca,'XTick',[-3:2]);
set(gca,'YTick',[-3:2]);
xlim([-3 2])
ylim([-3 2])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
plot([-3 2],[-3 2],'r')
xlabel(['Wake']);
ylabel(['REM']);

subplot(3,3,2)
hold on
plot(log10(StatesFR(isStruc&isccgUndef,swsfr)),log10(StatesFR(isStruc&isccgUndef, remfr)),'.','color',rgb('DarkGray'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgPyr,swsfr)),log10(StatesFR(isStruc&isccgPyr,remfr)),'.','color',rgb('OrangeRed'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgInt,swsfr)),log10(StatesFR(isStruc&isccgInt,remfr)),'.','color',rgb('MediumBlue'),'MarkerSize',15);
set(gca,'XTick',[-3:2]);
set(gca,'YTick',[-3:2]);
xlim([-3 2])
ylim([-3 2])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
plot([-3 2],[-3 2],'r')
xlabel(['SWS']);
ylabel(['REM']);

subplot(3,3,3)
hold on
plot(log10(StatesFR(isStruc&isccgUndef,wakefr)),log10(StatesFR(isStruc&isccgUndef,swsfr)),'.','color',rgb('DarkGray'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgPyr,swsfr)),log10(StatesFR(isStruc&isccgPyr,remfr)),'.','color',rgb('OrangeRed'),'MarkerSize',15);
plot(log10(StatesFR(isStruc&isccgInt,swsfr)),log10(StatesFR(isStruc&isccgInt,remfr)),'.','color',rgb('MediumBlue'),'MarkerSize',15);
set(gca,'XTick',[-3:2]);
set(gca,'YTick',[-3:2]);
xlim([-3 2])
ylim([-3 2])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
plot([-3 2],[-3 2],'r')
xlabel(['Wake']);
ylabel(['SWS']);

set (gcf,'Position',position);

meanFRccgint=mean(StatesFR(isStruc&isccgInt,[remfr swsfr wakefr]));
meanFRccgpyr=mean(StatesFR(isStruc&isccgPyr,[remfr swsfr wakefr]));
meanFRfinalint=mean(StatesFR(isStruc&isfinalInt,[remfr swsfr wakefr]));
meanFRfinalpyr=mean(StatesFR(isStruc&isfinalPyr,[remfr swsfr wakefr]));

semFRccgint=sem(StatesFR(isStruc&isccgInt,[remfr swsfr wakefr]));
semFRccgpyr=sem(StatesFR(isStruc&isccgPyr,[remfr swsfr wakefr]));
semFRfinalint=sem(StatesFR(isStruc&isfinalInt,[remfr swsfr wakefr]));
semFRfinalpyr=sem(StatesFR(isStruc&isfinalPyr,[remfr swsfr wakefr]));

subplot(3,3,6);
barwitherr([semFRfinalpyr 0 semFRccgpyr],[meanFRfinalpyr 0 meanFRccgpyr],'FaceColor',rgb('OrangeRed'),'LineStyle','none');
set(gca,'XTickLabel',{'REM' 'SWS' 'Wake' '/' 'REM' 'SWS' 'Wake'});
subplot(3,3,9);
barwitherr([semFRfinalint 0 semFRccgint],[meanFRfinalint 0 meanFRccgint],'FaceColor',rgb('MediumBlue'),'LineStyle','none');
set(gca,'XTickLabel',{'REM' 'SWS' 'Wake' '/' 'REM' 'SWS' 'Wake'});
xlabel('Mean state FR finalType/ccg-identified');

% gain distrib
allgains=StatesFR(:,remfr)./StatesFR(:,wakefr);
allgains(allgains==0)=NaN;

Mean.gain=nanmean(allgains(isStruc&isfinalPyr));
Sem.gain=nansem(allgains(isStruc&isfinalPyr));

allratios.rem=(StatesFR(:,remfr)-StatesFR(:,wakefr))./(StatesFR(:,remfr)+StatesFR(:,wakefr));
Mean.ratio=nanmean(allratios.rem(isStruc&isfinalPyr));
Sem.ratio=nansem(allratios.rem(isStruc&isfinalPyr));

allratios.sws=(StatesFR(:,swsfr)-StatesFR(:,wakefr))./(StatesFR(:,swsfr)+StatesFR(:,wakefr));

% Stats
skewn.pyr.all=skewness(log(allgains(isStruc&isfinalPyr)));
skewn.pyr.ccg=skewness(log(allgains(isStruc&isccgPyr)));
skewn.int.all=skewness(log(allgains(isStruc&isfinalInt)));
skewn.int.ccg=skewness(log(allgains(isStruc&isccgInt)));

%  psignrank.pyr.all=signrank(log(allgains(isStruc&isfinalPyr)));
%  psignrank.pyr.ccg=signrank(log(allgains(isStruc&isccgPyr)));
%  psignrank.int.all=signrank(log(allgains(isStruc&isfinalInt)));
%  psignrank.int.ccg=signrank(log(allgains(isStruc&isccgInt)));

psignrank.pyr.all=signrank(allratios.rem(isStruc&isfinalPyr));
[psignrank.pyr.ccg,h,stats]=signrank(allratios.rem(isStruc&isccgPyr))
psignrank.int.all=signrank(allratios.rem(isStruc&isfinalInt));
[psignrank.int.ccg,h,stats]=signrank(allratios.rem(isStruc&isccgInt))

pranksum.all=ranksum(allratios.rem(isStruc&isfinalPyr),allratios.rem(isStruc&isfinalInt));
pranksum.ccg=ranksum(allratios.rem(isStruc&isccgPyr),allratios.rem(isStruc&isccgInt));

[h1,b1]=hist(allratios.rem(isStruc&isfinalPyr),50);
[h2,b2]=hist(allratios.rem(isStruc&isfinalInt),b1);
[h3,b3]=hist(allratios.rem(isStruc&isccgPyr),b1);
[h4,b4]=hist(allratios.rem(isStruc&isccgInt),b1);
figure;hold on;
plot(b1,ZeroToOne(Smooth(h1,1)),'Color',rgb('OrangeRed'));
plot(b1,ZeroToOne(Smooth(h3,1)),'Color',rgb('Orange'));
plot(b1,ZeroToOne(Smooth(h2,1)),'Color',rgb('MediumBlue'));
plot(b1,ZeroToOne(Smooth(h4,1)),'Color',rgb('LightBlue'));
PlotHVLines(nanmedian(allratios.rem(isStruc&isfinalPyr)),'v','color',rgb('OrangeRed'));
PlotHVLines(nanmedian(allratios.rem(isStruc&isfinalInt)),'v','color',rgb('MediumBlue'));
PlotHVLines(nanmedian(allratios.rem(isStruc&isccgPyr)),'v','color',rgb('Orange'));
PlotHVLines(nanmedian(allratios.rem(isStruc&isccgInt)),'v','color',rgb('LightBlue'));
xlabel('REM/WAke')

[h1,b1]=hist(allratios.sws(isStruc&isfinalPyr),50);
[h2,b2]=hist(allratios.sws(isStruc&isfinalInt),b1);
[h3,b3]=hist(allratios.sws(isStruc&isccgPyr),b1);
[h4,b4]=hist(allratios.sws(isStruc&isccgInt),b1);
figure;hold on;
plot(b1,ZeroToOne(Smooth(h1,1)),'Color',rgb('OrangeRed'));
plot(b1,ZeroToOne(Smooth(h3,1)),'Color',rgb('Orange'));
plot(b1,ZeroToOne(Smooth(h2,1)),'Color',rgb('MediumBlue'));
plot(b1,ZeroToOne(Smooth(h4,1)),'Color',rgb('LightBlue'));
PlotHVLines(nanmedian(allratios.sws(isStruc&isfinalPyr)),'v','color',rgb('OrangeRed'));
PlotHVLines(nanmedian(allratios.sws(isStruc&isfinalInt)),'v','color',rgb('MediumBlue'));
PlotHVLines(nanmedian(allratios.sws(isStruc&isccgPyr)),'v','color',rgb('Orange'));
PlotHVLines(nanmedian(allratios.sws(isStruc&isccgInt)),'v','color',rgb('LightBlue'));
xlabel('SWS/WAke')


%  
med.pyr.all=nanmedian(log(allgains(isStruc&isfinalPyr)));
med.pyr.ccg=nanmedian(log(allgains(isStruc&isccgPyr)));
med.int.all=nanmedian(log(allgains(isStruc&isfinalInt)));
med.int.ccg=nanmedian(log(allgains(isStruc&isccgInt)));


figure(f0);
subplot(3,3,4:5);
[h1,bins]=hist(log(allgains(isStruc&isfinalPyr)),[-5:0.1:5]);
bar(bins,h1,'FaceColor',rgb('OrangeRed'),'LineStyle','none');
PlotHVLines(med.pyr.all,'v','k');
ylabel(['skewness ' num2str(skewn.pyr.all) ' p=' num2str(psignrank.pyr.all)]);
subplot(3,3,7:8);
h2 = hist(log(allgains(isStruc&isfinalInt)),bins);
bar(bins,h2,'FaceColor',rgb('MediumBlue'),'LineStyle','none');
PlotHVLines(med.int.all,'v','k');
xlabel('REM/Wake FR gain distribution - log scale');
ylabel(['skewness ' num2str(skewn.int.all) ' p=' num2str(psignrank.int.all)]);
suptitle(['All rats - State Firing Rates - ' strucname]);


