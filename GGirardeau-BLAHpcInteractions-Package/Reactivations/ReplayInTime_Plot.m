function ReplayInTime_Plot (struc,binsize,zsc,window,ctype,varargin)

%ReplayInTime_Plot - PLots Peri-ripple reactivation strength.
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
%    See also : ReplayInTime, ReplayInTime_All, ReplayInTime_LapTypes
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
  if strcmp(ctype,'all')
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat']);
  elseif strcmp(ctype,'pyr')
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat']);
  else	
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '.mat']);
  end 
elseif strcmp(stattype,'median')
  if strcmp(ctype,'all')
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '-MEDIAN.mat']);
  elseif strcmp(ctype,'pyr')
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '-MEDIAN.mat']);
  else	
    load(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-MEDIAN.mat']);
  end 
end  
  
  
f0=figure;
bar([mean(centermeans.pre.cross) mean(centermeans.post.cross)]);
hold on;
for i=1:length(centermeans.pre.cross)
  plot([1;2],[centermeans.pre.cross(i);centermeans.post.cross(i)]);
end

f1=figure('Position',[796 433 1123 339]);
subplot(1,3,1);
plot(tb,Mean.pre.cross,'b');hold on;
plot(tb,Mean.pre.cross+Sem.pre.cross,'b:');
plot(tb,Mean.pre.cross-Sem.pre.cross,'b:');
plot(tb,Mean.post.cross,'r');
plot(tb,Mean.post.cross+Sem.post.cross,'r:');
plot(tb,Mean.post.cross-Sem.post.cross,'r:');
xlabel(['Cross structure-' struc]);
ylabel('Reactivation strength');
ylim([-0.01 0.1])
subplot(1,3,2)
plot(tb,Mean.pre.hpc,'b');hold on;
plot(tb,Mean.post.hpc,'r');
plot(tb,Mean.pre.hpc+Sem.pre.hpc,'b:');hold on;
plot(tb,Mean.post.hpc+Sem.post.hpc,'r:');
plot(tb,Mean.pre.hpc-Sem.pre.hpc,'b:');hold on;
plot(tb,Mean.post.hpc-Sem.post.hpc,'r:');
xlabel('Intra-Hpc')
subplot(1,3,3)
plot(tb,Mean.pre.struc,'b');hold on;
plot(tb,Mean.post.struc,'r');
plot(tb,Mean.pre.struc+Sem.pre.struc,'b:');hold on;
plot(tb,Mean.post.struc+Sem.post.struc,'r:');
plot(tb,Mean.pre.struc-Sem.pre.struc,'b:');hold on;
plot(tb,Mean.post.struc-Sem.post.struc,'r:');
xlabel(['Intra-' struc]);
if strcmp(zsc,'on')
  suptitle(['Mean Raw Matrix z-scored reactivation strength n=' num2str(length(ratsess)) 'sess']);
else
  suptitle(['Mean Raw Matrix non z-scored reactivation strength' num2str(length(ratsess)) 'sess']);
end

f2=figure;
ratiodiff=Ratios.Post-Ratios.Pre;
meandiff=centermeans.post.cross-centermeans.pre.cross;
meandiff(isnan(ratiodiff))=[];
ratiodiff(isnan(ratiodiff))=[];
plot(ratiodiff,meandiff,'k.');
lsline
[coef,p]=corrcoef(ratiodiff,meandiff)

%  f2=figure('Position',[796 433 1123 339]);
%  subplot(1,3,1)
%  plot(tb,Cross.rat8.pre,'b');hold on;
%  plot(tb,Cross.rat8.post,'r');
%  plot(tb,Cross.rat8.pre+Cross.rat8.presem,'b:');
%  plot(tb,Cross.rat8.post+Cross.rat8.postsem,'r:');
%  plot(tb,Cross.rat8.pre-Cross.rat8.presem,'b:');
%  plot(tb,Cross.rat8.post-Cross.rat8.postsem,'r:');
%  xlabel('Rat8');
%  ylabel('Reactivation strength');
%  subplot(1,3,2)
%  plot(tb,Cross.rat10.pre,'b');hold on;
%  plot(tb,Cross.rat10.post,'r');
%  plot(tb,Cross.rat10.pre+Cross.rat10.presem,'b:');
%  plot(tb,Cross.rat10.post+Cross.rat10.postsem,'r:');
%  plot(tb,Cross.rat10.pre-Cross.rat10.presem,'b:');
%  plot(tb,Cross.rat10.post-Cross.rat10.postsem,'r:');
%  xlabel('Rat10');
%  subplot(1,3,3)
%  plot(tb,Cross.rat11.pre,'b');hold on;
%  plot(tb,Cross.rat11.post,'r');
%  plot(tb,Cross.rat11.pre+Cross.rat11.presem,'b:');
%  plot(tb,Cross.rat11.post+Cross.rat11.postsem,'r:');
%  plot(tb,Cross.rat11.pre-Cross.rat11.presem,'b:');
%  plot(tb,Cross.rat11.post-Cross.rat11.postsem,'r:');
%  xlabel('Rat11');
%  if strcmp(zsc,'on')
%    suptitle('Mean Raw Matrix z-scored reactivation strength-Per rat-Cross structure');
%  else
%    suptitle('Mean Raw Matrix non z-scored reactivation strength-PerRat-CrossStructure');
%  end
wind=tb<0.25&tb>-0.25;
wind=wind';

centermeans.pre.cross=mean(PeriRippleReplay.pre.cross(:,wind),2);
centermeans.post.cross=mean(PeriRippleReplay.post.cross(:,wind),2);
[h.cross,p.cross,stats]=signrank(centermeans.pre.cross,centermeans.post.cross,'tail','left')

keyboard

cd('/media/Data-01/All-Rats/AllRats-ReplayInTime/')
if strcmp(stattype,'median')
  plot2svg(['AllRats-PeriRippleReact-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '-MEDIAN.svg'],f1);
  plot2svg(['AllRats-PrePostMeanReact-CrossStruc-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '-MEDIAN.svg'],f0);
else
  plot2svg(['AllRats-PeriRippleReact-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '.svg'],f1);
  plot2svg(['AllRats-PrePostMeanReact-CrossStruc-' struc '-zsc' zsc '-binsize' num2str(binsize) '-window' int2str(window) '-ctype-' ctype '.svg'],f0);
end
%  plot2svg(['AllRats-PrePostMeanReact-CrossStruc-PerRat-' struc '-zsc' zsc '-binsize' num2str(binsize) '.svg'],f2);