function [PeriRippleReplay,Mean,Sem,tb,centermeans,ratsess] = ReplayInTime_LapTypes_All(struc,binsize,zsc,window,ctype,varargin)

%ReplayInTime_LapTypes_All - Reactivation strength of the safe vs airpuff RUNS in time during pre and post-sleep epochs across rats and sessions
%
%  USAGE
%
%   [PeriRippleReplay,Ratios,Mean,Sem,tb,centermeans,ratsess] = ReplayInTime_All(struc,binsize,zsc,window,ctype,varargin)
%
%    struc              structure name (ex : 'BLA')
%    binsize            binsize
%    zsc                zscore 'on'/'off'
%    window             periripple window
%    ctype              cell type 'pyr', 'all'
%    <options>
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'stattype'	'mean'(default), 'median'
%    =========================================================================
%
%  OUTPUT
%
%      PeriRippleReplay     peri-ripple reactivation strength (Matrix - one line/sess)
%      Ratios               BehavioralMeasures ratios to check for (non-existent) correlations 
%      Mean                 Mean prei-ripple replay across animals and sessions
%      Sem                  Sem for periripple replay acroos animals and sessions
%      tb                   timebins for peri-rippl replay
%      centermeans          Mean reactivation strength at the peak ripple (500ms window centered on ripple peak)
%      ratsess              lis of rat/sessions corresponding to the PeriRippleReplay matric and centermeans vector.
%
%  NOTE
%
%  SEE
%
%    See also : ReplayInTime, ReplayInTime_Plot, ReplayInTime_LapTypes
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

load('/media/Data-01/All-Rats/sessionindexing.mat');

PeriRippleReplay.pre.safe=[];
PeriRippleReplay.pre.ap=[];
PeriRippleReplay.post.safe=[];
PeriRippleReplay.post.ap=[];

ratsess=[];

for i=1:size(ratsessionindex)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  if exist ([xml '-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat'])==2;
    load([xml '-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat']);
    load([xml '-ExplainedVariance-' struc '.mat'],'CorrMatrix');
    if size(corrM.hpbla.safe,1)>=6 & size(corrM.hpbla.safe,2)>=6 
      ratsess=[ratsess;ratsessionindex(i,:)];
      PeriRippleReplay.pre.safe=[PeriRippleReplay.pre.safe;meanR.cross.pre.safe'];
      PeriRippleReplay.pre.ap=[PeriRippleReplay.pre.ap;meanR.cross.pre.ap'];    
      PeriRippleReplay.post.safe=[PeriRippleReplay.post.safe;meanR.cross.post.safe'];
      PeriRippleReplay.post.ap=[PeriRippleReplay.post.ap;meanR.cross.post.ap'];
    end
  end
end

if strcmp(stattype,'mean')
  Mean.pre.safe=mean(PeriRippleReplay.pre.safe,1);
  Mean.post.safe=mean(PeriRippleReplay.post.safe,1);
  Mean.pre.ap=mean(PeriRippleReplay.pre.ap,1);
  Mean.post.ap=mean(PeriRippleReplay.post.ap,1);
  Sem.pre.safe=sem(PeriRippleReplay.pre.safe);
  Sem.post.safe=sem(PeriRippleReplay.post.safe);
  Sem.pre.ap=sem(PeriRippleReplay.pre.ap);
  Sem.post.ap=sem(PeriRippleReplay.post.ap);
  
  MeanDiffTime.safe=mean(PeriRippleReplay.post.safe-PeriRippleReplay.pre.safe);
  MeanDiffTime.ap=mean(PeriRippleReplay.post.ap-PeriRippleReplay.pre.ap);
  SemDiffTime.safe=sem(PeriRippleReplay.post.safe-PeriRippleReplay.pre.safe);
  SemDiffTime.ap=sem(PeriRippleReplay.post.ap-PeriRippleReplay.pre.ap);
  
elseif strcmp(stattype,'median')
  Mean.pre.safe=median(PeriRippleReplay.pre.safe,1);
  Mean.post.safe=median(PeriRippleReplay.post.safe,1);
  Mean.pre.ap=median(PeriRippleReplay.pre.ap,1);
  Mean.post.ap=median(PeriRippleReplay.post.ap,1);
  Sem.pre.safe=semedian(PeriRippleReplay.pre.safe);
  Sem.post.safe=semedian(PeriRippleReplay.post.safe);
  Sem.pre.ap=semedian(PeriRippleReplay.pre.ap);
  Sem.post.ap=semedian(PeriRippleReplay.post.ap);
end  

wind=tb<0.25&tb>-0.25;
wind=wind';

centermeans.pre.safe=mean(PeriRippleReplay.pre.safe(:,wind),2);
centermeans.post.safe=mean(PeriRippleReplay.post.safe(:,wind),2);
[h.safe,p.safe]=signrank(centermeans.pre.safe,centermeans.post.safe,'tail','left');

centermeans.pre.ap=mean(PeriRippleReplay.pre.ap(:,wind),2);
centermeans.post.ap=mean(PeriRippleReplay.post.ap(:,wind),2);
[h.ap,p.ap]=signrank(centermeans.pre.ap,centermeans.post.ap,'tail','left')

nsess=length(ratsess)

cd('/media/Data-01/All-Rats/AllRats-ReplayInTime/');
if strcmp(stattype,'mean')
  save(['AllRats-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat'],'PeriRippleReplay','Mean','Sem','tb','centermeans','ratsess','MeanDiffTime','SemDiffTime');
elseif strcmp(stattype,'median')
  save(['AllRats-ReplayInTime-LapTypes-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '-MEDIAN.mat'],'PeriRippleReplay','Mean','Sem','tb','centermeans','ratsess');
end
