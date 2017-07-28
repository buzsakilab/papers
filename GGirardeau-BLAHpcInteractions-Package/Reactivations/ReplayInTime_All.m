function [PeriRippleReplay,Ratios,Mean,Sem,tb,centermeans,ratsess] = ReplayInTime_All(struc,binsize,zsc,window,ctype,varargin)

%ReplayInTime_All - Reactivation strength in time during pre and post-sleep epochs across rats and sessions
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

PeriRippleReplay.pre.cross=[];
PeriRippleReplay.post.cross=[];
PeriRippleReplay.pre.hpc=[];
PeriRippleReplay.post.hpc=[];
PeriRippleReplay.pre.struc=[];
PeriRippleReplay.post.struc=[];
Ratios.Pre=[];
Ratios.Post=[];

ratsess=[];

for i=1:size(ratsessionindex)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  if exist ([xml '-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat'])==2;
    load([xml '-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat']);
    load([xml '-ExplainedVariance-' struc '.mat'],'CorrMatrix');
    if exist ('BehavioralMeasures.mat')==2;
      load('BehavioralMeasures.mat');
    else
      ratio.prerun=NaN;
      ratio.postrun=NaN;
    end
    if strcmp(ctype,'pyr')
      if size(CorrMatrix.run.pyr.Rall,1)>=6 & size(CorrMatrix.run.pyr.Rall,2)>=6 
	ratsess=[ratsess;ratsessionindex(i,:)];
	PeriRippleReplay.pre.cross=[PeriRippleReplay.pre.cross;meanR.cross.pre'];
	PeriRippleReplay.post.cross=[PeriRippleReplay.post.cross;meanR.cross.post'];
	PeriRippleReplay.pre.hpc=[PeriRippleReplay.pre.hpc;meanR.hpc.pre'];
	PeriRippleReplay.post.hpc=[PeriRippleReplay.post.hpc;meanR.hpc.post'];
	PeriRippleReplay.pre.struc=[PeriRippleReplay.pre.struc;meanR.bla.pre'];
	PeriRippleReplay.post.struc=[PeriRippleReplay.post.struc;meanR.bla.post'];
	
	Ratios.Pre=[Ratios.Pre;ratio.prerun];
	Ratios.Post=[Ratios.Post;ratio.postrun];
  %        figure;
  %        subplot(1,2,1);hold on;
  %        plot(tb,meanR.cross.pre','b');
  %        plot(tb,meanR.cross.post','r');
  %        if exist('ratio')
  %  		subplot(1,2,2)
  %  		bar([ratio.prerun ratio.run ratio.postrun]);
  %        end
      end
    elseif strcmp(ctype,'all')
      if size(CorrMatrix.run.all.Rall,1)>=6 & size(CorrMatrix.run.all.Rall,2)>=6 
	ratsess=[ratsess;ratsessionindex(i,:)];
	PeriRippleReplay.pre.cross=[PeriRippleReplay.pre.cross;meanR.cross.pre'];
	PeriRippleReplay.post.cross=[PeriRippleReplay.post.cross;meanR.cross.post'];
	PeriRippleReplay.pre.hpc=[PeriRippleReplay.pre.hpc;meanR.hpc.pre'];
	PeriRippleReplay.post.hpc=[PeriRippleReplay.post.hpc;meanR.hpc.post'];
	PeriRippleReplay.pre.struc=[PeriRippleReplay.pre.struc;meanR.bla.pre'];
	PeriRippleReplay.post.struc=[PeriRippleReplay.post.struc;meanR.bla.post'];
	
	Ratios.Pre=[Ratios.Pre;ratio.prerun];
	Ratios.Post=[Ratios.Post;ratio.postrun];
  %        figure;
  %        subplot(1,2,1);hold on;
  %        plot(tb,meanR.cross.pre','b');
  %        plot(tb,meanR.cross.post','r');
  %        if exist('ratio')
  %  		subplot(1,2,2)
  %  		bar([ratio.prerun ratio.run ratio.postrun]);
  %        end
      end
    end
  end
end


figure;
subplot(1,2,1);hold on;
for i=1:size(PeriRippleReplay.pre.cross,1)
    plot(PeriRippleReplay.pre.cross(i,:))
end
subplot(1,2,2);hold on;
for i=1:size(PeriRippleReplay.post.cross,1)
    plot(PeriRippleReplay.post.cross(i,:))
end



if strcmp(stattype,'mean')
  Mean.pre.cross=mean(PeriRippleReplay.pre.cross,1);
  Mean.post.cross=mean(PeriRippleReplay.post.cross,1);
  Mean.pre.hpc=mean(PeriRippleReplay.pre.hpc,1);
  Mean.post.hpc=mean(PeriRippleReplay.post.hpc,1);
  Mean.pre.struc=mean(PeriRippleReplay.pre.struc,1);
  Mean.post.struc=mean(PeriRippleReplay.post.struc,1);
  Sem.pre.cross=sem(PeriRippleReplay.pre.cross);
  Sem.post.cross=sem(PeriRippleReplay.post.cross);
  Sem.pre.hpc=sem(PeriRippleReplay.pre.hpc);
  Sem.post.hpc=sem(PeriRippleReplay.post.hpc);
  Sem.pre.struc=sem(PeriRippleReplay.pre.struc);
  Sem.post.struc=sem(PeriRippleReplay.post.struc);
elseif strcmp(stattype,'median')
  Mean.pre.cross=median(PeriRippleReplay.pre.cross,1);
  Mean.post.cross=median(PeriRippleReplay.post.cross,1);
  Mean.pre.hpc=median(PeriRippleReplay.pre.hpc,1);
  Mean.post.hpc=median(PeriRippleReplay.post.hpc,1);
  Mean.pre.struc=median(PeriRippleReplay.pre.struc,1);
  Mean.post.struc=median(PeriRippleReplay.post.struc,1);
  Sem.pre.cross=semedian(PeriRippleReplay.pre.cross);
  Sem.post.cross=semedian(PeriRippleReplay.post.cross);
  Sem.pre.hpc=semedian(PeriRippleReplay.pre.hpc);
  Sem.post.hpc=semedian(PeriRippleReplay.post.hpc);
  Sem.pre.struc=semedian(PeriRippleReplay.pre.struc);
  Sem.post.struc=semedian(PeriRippleReplay.post.struc);
end 
 
 
wind=tb<0.25&tb>-0.25;
wind=wind';

centermeans.pre.cross=mean(PeriRippleReplay.pre.cross(:,wind),2);
centermeans.post.cross=mean(PeriRippleReplay.post.cross(:,wind),2);
[h.cross,p.cross]=signrank(centermeans.pre.cross,centermeans.post.cross,'tail','left');

centermeans.pre.hpc=mean(PeriRippleReplay.pre.hpc(:,wind),2);
centermeans.post.hpc=mean(PeriRippleReplay.post.hpc(:,wind),2);
[h.hpc,p.hpc]=signrank(centermeans.pre.hpc,centermeans.post.hpc,'tail','left')

centermeans.pre.struc=mean(PeriRippleReplay.pre.struc(:,wind),2);
centermeans.post.struc=mean(PeriRippleReplay.post.struc(:,wind),2);
[h.struc,p.struc]=signrank(centermeans.pre.struc,centermeans.post.struc,'tail','left');

nsess=length(ratsess)

%  Cross.rat8.pre=mean(PeriRippleReplay.pre.cross(ratsess(:,1)==8,:),1);
%  Cross.rat8.post=mean(PeriRippleReplay.post.cross(ratsess(:,1)==8,:),1);
%  Cross.rat10.pre=mean(PeriRippleReplay.pre.cross(ratsess(:,1)==10,:),1);
%  Cross.rat10.post=mean(PeriRippleReplay.post.cross(ratsess(:,1)==10,:),1);
%  Cross.rat11.pre=mean(PeriRippleReplay.pre.cross(ratsess(:,1)==11,:),1);
%  Cross.rat11.post=mean(PeriRippleReplay.post.cross(ratsess(:,1)==11,:),1);
%  Cross.rat8.presem=sem(PeriRippleReplay.pre.cross(ratsess(:,1)==8,:));
%  Cross.rat8.postsem=sem(PeriRippleReplay.post.cross(ratsess(:,1)==8,:));
%  Cross.rat10.presem=sem(PeriRippleReplay.pre.cross(ratsess(:,1)==10,:));
%  Cross.rat10.postsem=sem(PeriRippleReplay.post.cross(ratsess(:,1)==10,:));
%  Cross.rat11.presem=sem(PeriRippleReplay.pre.cross(ratsess(:,1)==11,:));
%  Cross.rat11.postsem=sem(PeriRippleReplay.post.cross(ratsess(:,1)==11,:));


%  cd('/media/Data-01/All-Rats/AllRats-ReplayInTime/');
%  if strcmp(stattype,'mean')
%    save(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '.mat'],'PeriRippleReplay','Ratios','Mean','Sem','tb','centermeans','ratsess');
%  elseif strcmp(stattype,'median')
%    save(['AllRats-ReplayInTime-' struc '-binsize' num2str(binsize) '-zsc' zsc '-window' int2str(window) '-ctype-' ctype '-MEDIAN.mat'],'PeriRippleReplay','Ratios','Mean','Sem','tb','centermeans','ratsess');
%  end
