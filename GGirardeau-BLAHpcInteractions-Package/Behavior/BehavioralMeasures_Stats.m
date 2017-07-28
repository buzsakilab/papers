function stats = BehavioralMeasures_Stats(varargin)

%  BehavioralMeasures_Stats - Statistical tests for behavioral measures.
%
%  USAGE
%
%   BehavioralMeasures_Stats
%
%  OUTPUT
%
%   Plots/P values
%
%  NOTE
%
%   Run BehavioralMeasures and BehavioralMeasuresAll beforehand
%
%  SEE
%  
%    BehavioralMeasures_Plot, BehavioralMeasures, BehavioralMeasuresAll
%
% Gabrielle Girardeau, March 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
ratnum='all';

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
    case 'ratnum',
	ratnum = varargin{i+1};
  otherwise,
	error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
  end
end

load('/media/Data-01/All-Rats/BehavioralMeasures/BehavioralMeasures.mat')

Rats=rats(:);
if strcmp(ratnum,'all')
  sesstotal=55;
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
else
  idx=find(Rats==ratnum)
  sesstotal=size(idx,1);
  for i=1:size(idx,1)
    prerun.currentvel(i)=PrePmeanVAll{idx(i)}.prerun;
    prerun.previousvel(i)=YPrePmeanVAll{idx(i)}.prerun;
    run.currentvel(i)=PrePmeanVAll{idx(i)}.run;
    run.previousvel(i)=YPrePmeanVAll{idx(i)}.run;
    postrun.currentvel(i)=PrePmeanVAll{idx(i)}.postrun;
    postrun.previousvel(i)=YPrePmeanVAll{idx(i)}.postrun;
    ratios.prerun(i)=ratioAll{idx(i)}.prerun;
    ratios.run(i)=ratioAll{idx(i)}.run;
    ratios.postrun(i)=ratioAll{idx(i)}.postrun;
  end
end


%% Test normality
[h,p]=kstest(ratios.prerun)
[h,p]=kstest(ratios.run)
[h,p]=kstest(ratios.run)
%  [h,p]=kstest(prerun.currentvel)
%  [h,p]=kstest(prerun.previousvel)
%  [h,p]=kstest(run.previousvel)
%  [h,p]=kstest(run.currentvel)
%  [h,p]=kstest(postrun.currentvel)
%  [h,p]=kstest(postrun.previousvel)

%% Test equality of variances
[a,b]=vartestn([ratios.prerun' ratios.run' ratios.postrun'])
%  [a,b]=vartestn([prerun.currentvel' prerun.previousvel' run.previousvel' run.currentvel' postrun.currentvel' postrun.previousvel'])

%%%%% Two-way anova on velocities : factor 1 PRE/RUN/POST, factor2 current/previous AP loc.
isNaN=isnan(prerun.currentvel)|isnan(prerun.previousvel)|isnan(run.currentvel)|isnan(run.previousvel)|isnan(postrun.currentvel)|isnan(postrun.previousvel);
a=sum(isNaN)
siz=sesstotal-a;
y=[prerun.currentvel(~isNaN)';prerun.previousvel(~isNaN)';run.currentvel(~isNaN)';run.previousvel(~isNaN)';postrun.currentvel(~isNaN)';postrun.previousvel(~isNaN)'];
subj=repmat([1:siz]',[6 1]);
f1=[ones(2*siz,1);ones(2*siz,1)*2;ones(2*siz,1)*3];
f3=[ones(siz,1);ones(siz,1)*2];
f2=[f3;f3;f3];
stats.vel=rm_anova2(y,subj,f1,f2,{'prerunpost', 'currprev'});

%%%%%% Manual post-hoc tests (paired ttests)
[h1,p1,~,stats1]=ttest(prerun.currentvel,prerun.previousvel);
[h2,p2,~,stats2]=ttest(run.currentvel,run.previousvel);
[h3,p3,~,stats3]=ttest(postrun.currentvel,postrun.previousvel);
[h4,p4]=ttest(prerun.currentvel,postrun.currentvel);

%%%%%% One way repeated measure anova on ratios + post-hoc comparison.
isNaN2=isnan(ratios.prerun)|isnan(ratios.run)|isnan(ratios.postrun);
b=sum(isNaN2);
siz2=sesstotal-b;
zz=[ratios.prerun(~isNaN2)' ratios.run(~isNaN2)' ratios.postrun(~isNaN2)'];
[p,table,stats.ratios] = anova_rm(zz)
%  [comp,m]=multcompare(stats.ratios,'CType','bonferroni','Alpha',0.05)
[comp,m]=multcompare(stats.ratios,'Alpha',0.001)

[h5,p5,~,st5]=ttest(ratios.prerun,ratios.run)
[h6,p6,~,st6]=ttest(ratios.prerun,ratios.postrun)
[h7,p7,~,st7]=ttest(ratios.run,ratios.postrun)

figure;
subplot(1,2,1)
boxplot([ratios.prerun' ratios.run' ratios.postrun'])
hold on;
plot(1,ratios.prerun,'k.');
plot(2,ratios.run,'k.');
plot(3,ratios.postrun,'k.');
plot([nanmean(ratios.prerun) nanmean(ratios.run) nanmean(ratios.postrun)],'g.','MarkerSize',10)
subplot(1,2,2)
boxplot([prerun.previousvel' prerun.currentvel' run.previousvel' run.currentvel' postrun.previousvel' postrun.currentvel'])
hold on;
plot(1,prerun.previousvel,'k.');
plot(2,prerun.currentvel,'k.');
plot(3,run.previousvel,'k.');
plot(4,run.currentvel,'k.');
plot(5,postrun.previousvel,'k.');
plot(6,postrun.currentvel,'k.');
plot([nanmean(prerun.previousvel) nanmean(prerun.currentvel) nanmean(run.previousvel) nanmean(run.currentvel) nanmean(postrun.previousvel) nanmean(postrun.currentvel)],'g.','MarkerSize',10)
xlabel(['Rat' int2str(ratnum)])
