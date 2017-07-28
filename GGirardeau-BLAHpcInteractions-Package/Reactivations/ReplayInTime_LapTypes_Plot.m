function ReplayInTime_LapTypes_Plot (struc,binsize,zsc,window,ctype,varargin)

%ReplayInTime_LapTypes_Plot - PLots Peri-ripple reactivation strength for airpuff and safe RUN trajectories.
%
%  USAGE
%
%    ReplayInTime_Plot (struc,binsize,zsc,window,ctype,varargin)
%
%    struc              structure name (ex : 'BLA')
%    binsize            binsize
%    zsc                zscore 'on'/'off'
%    window             periripple window
%    ctype              cell type 'pyr', 'all'
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'stattype'	'mean' (default) or 'median'
%    =========================================================================
%
%  OUTPUT
%
%    PLots
%
%  NOTE
%  
%   INput variables are used to select pre-calculated variable, that must have been obtained beforehand with the same parameters using ReplayInTime.m/ReplayInTime_All
%
%  SEE
%
%    See also : ReplayInTime_LapTypes, ReplayInTime_LapTypes_All, ReplayInTime
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
stattype = 'mean';

% Check number of inputs
if nargin < 5,
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
    case 'stattype',
      stattype = lower(varargin{i+1});
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
   end
end

cd('/media/Data-01/All-Rats/AllRats-ReplayInTime/');
if strcmp(stattype,'mean')
  load(['AllRats-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat']);
elseif strcmp(stattype,'median')
  load(['AllRats-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '-MEDIAN.mat']);
end 


f0=figure;
subplot(1,2,1)
%  bar([mean(centermeans.pre.safe) mean(centermeans.post.safe)]);
boxplot([centermeans.pre.safe centermeans.post.safe],'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2);
hold on;
for i=1:length(centermeans.pre.safe)
  plot([1;2],[centermeans.pre.safe(i);centermeans.post.safe(i)]);
end
xlabel('safe')
subplot(1,2,2)
%  bar([mean(centermeans.pre.ap) mean(centermeans.post.ap)]);
boxplot([centermeans.pre.ap centermeans.post.ap],'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2);
ylim([-0.03 0.17]);
hold on;
for i=1:length(centermeans.pre.ap)
  plot([1;2],[centermeans.pre.ap(i);centermeans.post.ap(i)]);
end
ylim([-0.03 0.17])
xlabel('airpuff')

f1=figure('Position',[796 433 1123 339]);
subplot(1,2,1);
plot(tb,Mean.pre.safe,'b');hold on;
plot(tb,Mean.pre.safe+Sem.pre.safe,'b:');
plot(tb,Mean.pre.safe-Sem.pre.safe,'b:');
plot(tb,Mean.post.safe,'r');
plot(tb,Mean.post.safe+Sem.post.safe,'r:');
plot(tb,Mean.post.safe-Sem.post.safe,'r:');
ylim([-0.01 0.05])
xlabel(['SAFE-' struc]);
ylabel('Reactivation strength');
subplot(1,2,2);
plot(tb,Mean.pre.ap,'b');hold on;
plot(tb,Mean.pre.ap+Sem.pre.ap,'b:');
plot(tb,Mean.pre.ap-Sem.pre.ap,'b:');
plot(tb,Mean.post.ap,'r');
plot(tb,Mean.post.ap+Sem.post.ap,'r:');
plot(tb,Mean.post.ap-Sem.post.ap,'r:');
xlabel(['AIRPUFF-' struc]);
ylabel('Reactivation strength');
ylim([-0.01 0.05])

if strcmp(zsc,'on')
  suptitle(['Mean Raw Matrix z-scored reactivation strength n=' num2str(length(ratsess)) 'sess']);
else
  suptitle(['Mean Raw Matrix non z-scored reactivation strength' num2str(length(ratsess)) 'sess']);
end


%  f0=figure('Position',[796 433 1123 339]);hold on;
%  plot(tb,MeanDiffTime.safe,'b');hold on;
%  plot(tb,MeanDiffTime.safe+SemDiffTime.safe,'b:');
%  plot(tb,MeanDiffTime.safe-SemDiffTime.safe,'b:');
%  ylim([-0.01 0.05]);
%  ylabel('Reactivation strength difference');
%  plot(tb,MeanDiffTime.ap,'r');hold on;
%  plot(tb,MeanDiffTime.ap+SemDiffTime.ap,'r:');
%  plot(tb,MeanDiffTime.ap-SemDiffTime.ap,'r:');

meandiff.safe=centermeans.post.safe-centermeans.pre.safe;
meandiff.ap=centermeans.post.ap-centermeans.pre.ap;
[hdiff,pdiff,stdiff]=signrank(meandiff.safe,meandiff.ap,'tail','left')

f2=figure;hold on;
boxplot([meandiff.safe meandiff.ap],'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2);
for i=1:length(meandiff.safe)
    plot([1 2],[meandiff.safe(i) meandiff.ap(i)],'k');
end
xlabel('Pre/post Rd peak diff')
set(gca,'XTickLabel',{'Safe' 'Airpuff'})

wind=tb<0.25&tb>-0.25;
wind=wind';

centermeans.pre.safe=mean(PeriRippleReplay.pre.safe(:,wind),2);
centermeans.pre.ap=mean(PeriRippleReplay.pre.ap(:,wind),2);
centermeans.post.safe=mean(PeriRippleReplay.post.safe(:,wind),2);
centermeans.post.ap=mean(PeriRippleReplay.post.ap(:,wind),2);

[h.safe,p.safe,stats.safe]=signrank(centermeans.pre.safe,centermeans.post.safe,'tail','left')
[h.ap,p.ap,stats.ap]=signrank(centermeans.pre.ap,centermeans.post.ap,'tail','left')

keyboard

%  cd('/media/Data-01/All-Rats/AllRats-ReplayInTime/')
%  if strcmp(stattype,'mean')
%    plot2svg(['AllRats-PeriRippleReact-LapTypes' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '.svg'],f1);
%    plot2svg(['AllRats-PrePostMeanReact-LapTypes-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '.svg'],f0);
%    plot2svg(['AllRats-RDiffBoxplot-LapTypes-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '.svg'],f2);
%  else
%    plot2svg(['AllRats-PeriRippleReact-LapTypes' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '-MEDIAN.svg'],f1);
%    plot2svg(['AllRats-PrePostMeanReact-LapTypes-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '-MEDIAN.svg'],f0);
%    plot2svg(['AllRats-RDiffBoxplot-LapTypes-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '-MEDIAN.svg'],f2);
%  end