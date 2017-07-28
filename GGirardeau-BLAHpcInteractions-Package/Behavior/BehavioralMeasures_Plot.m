function BehavioralMeasures_Plot

%  BehavioralMeasures_Plot -PLots Mean speed in the current/previous danger zone for preRUN, RUN and postRUN epochs and associated Ratios across rats and sessions
%
%  USAGE
%
%   BehavioralMeasures_Plot 
%
%  OUTPUT
%
%   Plots
%
%  NOTE
%
%   Run BehavioralMeasures and BehavioralMeasuresAll beforehand
%
%  SEE
%  
%    BehavioralMeasures_Stats, BehavioralMeasures, BehavioralMeasuresAll
%
% Gabrielle Girardeau, March 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/BehavioralMeasures/BehavioralMeasures.mat')
load(['/media/Data-01/All-Rats/AllRats-ReplayInTime/AllRats-ReplayInTime-LapTypes-BLA-binsize0.05-zscon-window2-ctype-pyr.mat']);

% hardcoded for number of sessions
for i=1:55
  prerun.currentvel(i)=PrePmeanVAll{i}.prerun;
  prerun.previousvel(i)=YPrePmeanVAll{i}.prerun;
  run.currentvel(i)=PrePmeanVAll{i}.run;
  run.previousvel(i)=YPrePmeanVAll{i}.run;
  postrun.currentvel(i)=PrePmeanVAll{i}.postrun;
  postrun.previousvel(i)=YPrePmeanVAll{i}.postrun;
  ratios.prerun(i)=ratioAll{i}.prerun;
  ratios.run(i)=ratioAll{i}.run;
  ratios.postrun(i)=ratioAll{i}.postrun;
end

a=ones(55,1);
b=a*2;
c=a*3;
d=a*4;
e=a*5;
f=a*6;

%  figure;
%  barwitherr([nansem(ratios.prerun) nansem(ratios.run) nansem(ratios.postrun)],[nanmean(ratios.prerun) nanmean(ratios.run) nanmean(ratios.postrun)]);
%  hold on
%  plot(a,ratios.prerun,'r.');
%  plot(b,ratios.run,'r.');
%  plot(c,ratios.postrun,'r.');

figure;
boxplot([ratios.prerun' ratios.run' ratios.postrun'],'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2)


%  figure;
%  barwitherr([nansem(prerun.currentvel) nansem(run.currentvel) nansem(postrun.currentvel)],[nanmean(prerun.currentvel) nanmean(run.currentvel) nanmean(postrun.currentvel)]);
%  [h,p]=ttest2(prerun.currentvel,postrun.currentvel)
%  hold on
%  plot(a,prerun.currentvel,'r.');
%  plot(b,run.currentvel,'r.');
%  plot(c,postrun.currentvel,'r.');
%  
%  figure;
%  barwitherr([nansem(prerun.previousvel) nansem(run.previousvel) nansem(postrun.previousvel)],[nanmean(prerun.previousvel) nanmean(run.previousvel) nanmean(postrun.previousvel)]);
%  hold on
%  plot(a,prerun.previousvel,'r.');
%  plot(b,run.previousvel,'r.');
%  plot(c,postrun.previousvel,'r.');

figure;
barwitherr([nansem(prerun.previousvel) nansem(prerun.currentvel) nansem(run.previousvel) nansem(run.currentvel) nansem(postrun.previousvel) nansem(postrun.currentvel)],[nanmean(prerun.previousvel) nanmean(prerun.currentvel) nanmean(run.previousvel) nanmean(run.currentvel) nanmean(postrun.previousvel) nanmean(postrun.currentvel)]);
hold on
plot(a,prerun.previousvel,'k.');
plot(c,run.previousvel,'k.');
plot(e,postrun.previousvel,'k.');
plot(b,prerun.currentvel,'r.');
plot(d,run.currentvel,'r.');
plot(f,postrun.currentvel,'r.');

figure;
boxplot([prerun.previousvel' prerun.currentvel' run.previousvel' run.currentvel' postrun.previousvel' postrun.currentvel'],'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2);

%  ratios.prerun
%  [ratios.prerun ratios.run ratios.postrun]

[p,tab,stat]=anova1([ratios.prerun' ratios.run' ratios.postrun'])
multcompare(stat,'alpha',0.001)

[h,pt]=ttest(ratios.prerun,ratios.postrun)
[h,pt]=ttest(prerun.currentvel,postrun.currentvel)
[h,pt]=ttest(postrun.currentvel,postrun.previousvel)

%%%%%%%%%%%%%%%%%%%%%%%%%%








