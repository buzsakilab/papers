function [r,p] = PrePostFiringRates_Plot(CellList,varargin)

%PrePostFiringRates_Plots - Plots firing rates before vs after learning for a definite group of cells
%
%  USAGE
%  
%    r = PrePostFiringRates_Plot(CellList,varargin)
%
%    CellList           List of cells to include in the plot [Rat / sess/ shank / unit / id]
%    <options>          List of property-value pairs
%  
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'newfig'		plot in new figure 'on'/'off'
%     'color'		color for plot ('k' default)
%     'savefig'         save figure ('on'/'off')
%     'stattest'        Statistical test (Ttest [default], SignRank)
%    =========================================================================
%
%  OUTPUT
%
%    r      correlation coefficient for pre/post rat correlation
%    p      pvalue for r
%
%  NOTES
%
%  SEE
%
% 2016 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
newfig= 'on';
color='k';
savefig='off';
stattest='Ttest'; % or 'SignRank'

% Check varargin
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters  (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
  switch(lower(varargin{i})),
    case 'newfig'
      newfig = varargin{i+1};
    case 'color'
      color= varargin{i+1};
    case 'savefig'
      savefig= varargin{i+1};
    case 'stattest'
      stattest = varargin{i+1};
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

load('/media/Data-01/All-Rats/AllRats-PrePostFRAll.mat');

prerem=6;
postrem=7;
presws=8;
postsws=9;
prerun=12;
run=13;
postrun=14;

CellFR=PrePostFRAll(ismember(PrePostFRAll(:,1:4),CellList,'rows'),:);

cellgroup = input('CellGroup?','s');

%%%%%%%%%%%% Stats
if strcmp(stattest,'SignRank')
  [prun,~]=signrank(CellFR(:,prerun),CellFR(:,postrun))
  [psws,~]=signrank(CellFR(:,presws),CellFR(:,postsws))
  [prem,~]=signrank(CellFR(:,prerem),CellFR(:,postrem))
elseif strcmp(stattest,'Ttest')
  [~,prun]=ttest(CellFR(:,prerun),CellFR(:,postrun))
  [~,psws]=ttest(CellFR(:,presws),CellFR(:,postsws))
  [~,prem]=ttest(CellFR(:,prerem),CellFR(:,postrem))
end

%%%%%%%%%%% Plot
if strcmp(newfig,'on')
  figure('Position',[882 314 895 585]);
end
subplot(2,3,1);hold on;
plot(log10(CellFR(:,prerun)),log10(CellFR(:,postrun)),'.','Color',rgb('Crimson'),'MarkerSize',15);
%  pf=polyfit(log10(CellFR(:,prerun)),log10(CellFR(:,postrun)),1);
%  yfit=polyval(pf,[-3.5 1.5]);
%  plot([-3.5 1.5],yfit,'Color',color);
plot([-3.5 1.5],[-3.5 1.5],'k:')
[r.run, p.run]=corrcoef(CellFR(:,prerun),CellFR(:,postrun));
ylim([-3.5 1.5])
xlim([-3.5 1.5])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
xlabel('preRUN FR');
ylabel('postRUN FR');
%  title(['r=' num2str(p.run(1,2)) ' p=' num2str(p.run(1,2))]);
title([stattest ' p=' num2str(prun)]);
subplot(2,3,4);hold on;
bar([-1 1],[mean(CellFR(:,prerun)) mean(CellFR(:,postrun))],'FaceColor',rgb('Crimson'))
xlim([-5 5]);

subplot(2,3,2)
hold on
plot(log10(CellFR(:,prerem)),log10(CellFR(:,postrem)),'.','Color',rgb('SeaGreen'),'MarkerSize',15);
plot([-3.5 1.5],[-3.5 1.5],'k:')
%  lsline(gca);
[r.rem, p.rem]=corrcoef(CellFR(:,prerem),CellFR(:,postrem));
ylim([-3.5 1.5])
xlim([-3.5 1.5])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
xlabel('preREM FR');
ylabel('postREM FR');
%  title(['r=' num2str(p.rem(1,2)) ' p=' num2str(p.rem(1,2))]);
title([stattest ' p=' num2str(prem)]);
subplot(2,3,5);hold on;
bar([-1 1],[mean(CellFR(:,prerem)) mean(CellFR(:,postrem))],'FaceColor',rgb('SeaGreen'))
xlim([-5 5]);

subplot(2,3,3)
hold on
plot([-3.5 1.5],[-3.5 1.5],'k:')
plot(log10(CellFR(:,presws)),log10(CellFR(:,postsws)),'.','Color',rgb('SteelBlue'),'MarkerSize',15);
ylim([-3.5 1.5])
xlim([-3.5 1.5])
set(gca,'XTickLabel',[0.001 0.01 0.1 1 10 100]);
set(gca,'YTickLabel',[0.001 0.01 0.1 1 10 100]);
%  lsline(gca);
[r.sws,p.sws]=corrcoef(CellFR(:,presws),CellFR(:,postsws));
xlabel('preSWS FR');
ylabel('postSWS FR');
%  title(['r=' num2str(p.sws(1,2)) ' p=' num2str(p.sws(1,2))]);
title([stattest ' p=' num2str(psws)]);
subplot(2,3,6);hold on;
bar([-1 1],[mean(CellFR(:,presws)) mean(CellFR(:,postsws))],'FaceColor',rgb('Steelblue'))
xlim([-5 5]);

suptitle(['Pre vs Post FR (log10)- ' cellgroup]);

if strcmp(savefig,'on')
  cd('/media/Data-01/All-Rats/AllRats-PrePostFR')
  plot2svg(['PrePostFR-' cellgroup '.svg'],gcf);
end
  