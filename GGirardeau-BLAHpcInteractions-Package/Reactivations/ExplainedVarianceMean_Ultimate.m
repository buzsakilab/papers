function [sess,EVREV,p,stats] = ExplainedVarianceMean_Ultimate(struct,varargin)

% ExplainedVarianceMean_Ultimate - Calculates the mean explained variance and reverse explained variances across sessions (Kudrimoti 1999, Lansink 2009)
% ExplainedVariance_Ultimate must have run first
%  
%  USAGE
%
%    [sess, EVREV] = ExplainedVarianceMean_Ultimate (runsessionnumber,airpuffXlocation,cell)
%
%    structure		structure name (ex : 'BLA')
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%    EVtype		'cross' (cross-structure : default), 'intra' (intrastructure) 			
%    =========================================================================
%
%  OUTPUT
%
%    EVREV           	structure of vectors [EV REV], for SWS [EV REV EV1 REV1 EV2 REV2 EV3 REV3] 
%    sess		session list
%
%  NOTE
%
%
%  SEE
%
%    See also : binspikes, corrcoef, corr
%
% December 2015, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
EVtype='cross';
skiprip=0;
savefig='off';
savevar='off';

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
	  case 'evtype',
		EVtype = lower(varargin{i+1});
	  case 'skiprip',
		skiprip = lower(varargin{i+1});
	  case 'savevar'
		savevar = lower(varargin{i+1});
	  otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end



load('/media/Data-01/All-Rats/sessionindexing.mat')

EVREV.sws.all.Rall=[];
EVREV.sws.all.Rmod=[];
EVREV.sws.all.Rnomod=[];

EVREV.rem.all.Rall=[];
EVREV.rem.all.Rmod=[];
EVREV.rem.all.Rnomod=[];

EVREV.rip_in.all.Rall=[];
EVREV.rip_out.all.Rall=[];
EVREV.rip_in.all.Rmod=[];
EVREV.rip_in.all.Rnomod=[];
EVREV.rip_out.all.Rmod=[];
EVREV.rip_out.all.Rnomod=[];

EVREV.sws.pyr.Rall=[];
EVREV.sws.pyr.Rmod=[];
EVREV.sws.pyr.Rnomod=[];

EVREV.rem.pyr.Rall=[];
EVREV.rip_in.pyr.Rall=[];
EVREV.rip_out.pyr.Rall=[];

EVREV.rem.pyr.Rmod=[];
EVREV.rem.pyr.Rnomod=[];
EVREV.rip_in.pyr.Rmod=[];
EVREV.rip_in.pyr.Rnomod=[];
EVREV.rip_out.pyr.Rmod=[];
EVREV.rip_out.pyr.Rnomod=[];

EVREV.run.all.Rall=[];
EVREV.run.all.Rmod=[];
EVREV.run.all.Rnomod=[];
EVREV.run.pyr.Rall=[];
EVREV.run.pyr.Rmod=[];
EVREV.run.pyr.Rnomod=[];

sess.sws.all.Rall=[];
sess.rem.all.Rall=[];
sess.rem.all.Rmod=[];
sess.sws.all.Rmod=[];
sess.rip.all.Rall=[];
sess.rip.all.Rmod=[];
sess.sws.pyr.Rall=[];
sess.rem.pyr.Rall=[];
sess.rem.pyr.Rmod=[];
sess.sws.pyr.Rmod=[];
sess.rip.pyr.Rall=[];
sess.rip.pyr.Rmod=[];

sess.run.all.Rall=[];
sess.run.all.Rmod=[];
sess.run.pyr.Rall=[];
sess.run.pyr.Rmod=[];


for i=1:length(xmlpath)
  xml=xmlpath{i}(end-14:end-1);
  cd(xmlpath{i})
  ratsess=ratsessionindex(strcmp(xmlpath{i},xmlpath),1:2); %[rat session]
  if exist('BehavioralMeasures.mat')==2
    load('BehavioralMeasures.mat')
  else
    ratio.run=NaN;
    ratio.postrun=NaN;
  end  
  if exist('Perf.mat')==2
    load('Perf.mat')
  else
    perf=NaN;
  end
  
  if strcmp(EVtype,'cross') & exist([xml '-ExplainedVariance-' struct '.mat'])==2
    load([xml '-ExplainedVariance-' struct '.mat'])
    treatsess=1;
    if isfield(CorrMatrix.pre,'run')
      treatrun=1;
    else
      treatrun=0;
    end
  elseif strcmp(EVtype,'intra') & exist([xml '-ExplainedVariance-' struct '-IntraStructure.mat'])==2
    load([xml '-ExplainedVariance-' struct '-IntraStructure.mat']);
    treatsess=1;
    treatrun=0;
  else
    treatsess=0;
    treatrun=0;
  end

  % Minimum 6 cells in each structure (or intrastruct)
  enough.all.Rall=size(CorrMatrix.pre.sws.all.Rall,1)>5&size(CorrMatrix.pre.sws.all.Rall,2)>5;
  if isfield(CorrMatrix.pre.sws.all,'Rmod')
    enough.all.Rmod=size(CorrMatrix.pre.sws.all.Rmod,1)>5&size(CorrMatrix.pre.sws.all.Rmod,2)>5;
  else
    enough.all.Rmod=0;
  end
  
  if isfield(CorrMatrix.pre.sws,'pyr')
    enough.pyr.Rall=size(CorrMatrix.pre.sws.pyr.Rall,1)>5&size(CorrMatrix.pre.sws.pyr.Rall,2)>5;
    if isfield(CorrMatrix.pre.sws.pyr,'Rmod')
      enough.pyr.Rmod=size(CorrMatrix.pre.sws.pyr.Rmod,1)>5&size(CorrMatrix.pre.sws.pyr.Rmod,2)>5;
    else
      enough.pyr.Rmod=0;
    end  
  else
    enough.pyr.Rall=0;
    enough.pyr.Rmod=0;
  end
    
  enough.rem=duration.pre.rem>180 & duration.post.rem>180;
  enough.rippledata=isfield(EV,'rip_in');
  
  if treatsess
    %% Alltypes, Allripplemod
    if enough.all.Rall
      EVREV.sws.all.Rall=[EVREV.sws.all.Rall; EV.sws.all.Rall REV.sws.all.Rall EV.post1.sws.all.Rall REV.post1.sws.all.Rall EV.post2.sws.all.Rall REV.post2.sws.all.Rall EV.post3.sws.all.Rall REV.post3.sws.all.Rall];    
      sess.sws.all.Rall=[sess.sws.all.Rall;xml];
      if enough.rippledata
	EVREV.rip_in.all.Rall=[EVREV.rip_in.all.Rall; EV.rip_in.all.Rall REV.rip_in.all.Rall];
	EVREV.rip_out.all.Rall=[EVREV.rip_out.all.Rall; EV.rip_out.all.Rall REV.rip_out.all.Rall];
	sess.rip.all.Rall=[sess.rip.all.Rall;xml];
      end
      if enough.rem
	EVREV.rem.all.Rall=[EVREV.rem.all.Rall; EV.rem.all.Rall REV.rem.all.Rall];
	sess.rem.all.Rall=[sess.rem.all.Rall;xml];
      end
      if treatrun
	EVREV.run.all.Rall=[EVREV.run.all.Rall;EV.run.all.Rall REV.run.all.Rall];
	sess.run.all.Rall=[sess.run.all.Rall;xml];
      end
    end
    %% Alltypes, Ripple Mod vs Nomod
    if enough.rippledata
      if enough.all.Rmod % Minimum 6 cells in each structure
	EVREV.sws.all.Rmod=[EVREV.sws.all.Rmod; EV.sws.all.Rmod REV.sws.all.Rmod EV.post1.sws.all.Rmod REV.post1.sws.all.Rmod EV.post2.sws.all.Rmod REV.post2.sws.all.Rmod EV.post3.sws.all.Rmod REV.post3.sws.all.Rmod];
	EVREV.sws.all.Rnomod=[EVREV.sws.all.Rnomod; EV.sws.all.Rnomod REV.sws.all.Rnomod EV.post1.sws.all.Rnomod REV.post1.sws.all.Rnomod EV.post2.sws.all.Rnomod REV.post2.sws.all.Rnomod EV.post3.sws.all.Rnomod REV.post3.sws.all.Rnomod];	
	EVREV.rip_in.all.Rmod=[EVREV.rip_in.all.Rmod; EV.rip_in.all.Rmod REV.rip_in.all.Rmod];
	EVREV.rip_in.all.Rnomod=[EVREV.rip_in.all.Rnomod; EV.rip_in.all.Rnomod REV.rip_in.all.Rnomod];
	EVREV.rip_out.all.Rmod=[EVREV.rip_out.all.Rmod; EV.rip_out.all.Rmod REV.rip_out.all.Rmod];
	EVREV.rip_out.all.Rnomod=[EVREV.rip_out.all.Rnomod; EV.rip_out.all.Rnomod REV.rip_out.all.Rnomod];
	sess.sws.all.Rmod=[sess.sws.all.Rmod;xml];
	sess.rip.all.Rmod=[sess.rip.all.Rmod;xml];
	if enough.rem
	  EVREV.rem.all.Rmod=[EVREV.rem.all.Rmod; EV.rem.all.Rmod REV.rem.all.Rmod];
	  EVREV.rem.all.Rnomod=[EVREV.rem.all.Rnomod; EV.rem.all.Rnomod REV.rem.all.Rnomod];
	  sess.rem.all.Rmod=[sess.rem.all.Rmod;xml];
	end
	if treatrun
	  EVREV.run.all.Rmod=[EVREV.run.all.Rmod; EV.run.all.Rmod REV.run.all.Rmod];
	  EVREV.run.all.Rnomod=[EVREV.run.all.Rnomod; EV.run.all.Rnomod REV.run.all.Rnomod];
	  sess.run.all.Rmod=[sess.run.all.Rmod;xml];
	end
      end
    end
    %% Pyr Only, Allripplemod
    if enough.pyr.Rall % Minimum 6 cells in each structure
      EVREV.sws.pyr.Rall=[EVREV.sws.pyr.Rall; EV.sws.pyr.Rall REV.sws.pyr.Rall EV.post1.sws.pyr.Rall REV.post1.sws.pyr.Rall EV.post2.sws.pyr.Rall REV.post2.sws.pyr.Rall EV.post3.sws.pyr.Rall REV.post3.sws.pyr.Rall];
      sess.sws.pyr.Rall=[sess.sws.pyr.Rall;xml];
      if enough.rem
	EVREV.rem.pyr.Rall=[EVREV.rem.pyr.Rall; EV.rem.pyr.Rall REV.rem.pyr.Rall];
	sess.rem.pyr.Rall=[sess.rem.pyr.Rall;xml];
      end
      if treatrun
	EVREV.run.pyr.Rall=[EVREV.run.pyr.Rall; EV.run.pyr.Rall REV.run.pyr.Rall];
	sess.run.pyr.Rall=[sess.run.pyr.Rall;xml];
      end
      if enough.rippledata
	EVREV.rip_in.pyr.Rall=[EVREV.rip_in.pyr.Rall; EV.rip_in.pyr.Rall REV.rip_in.pyr.Rall];
	EVREV.rip_out.pyr.Rall=[EVREV.rip_out.pyr.Rall; EV.rip_out.pyr.Rall REV.rip_out.pyr.Rall];
	sess.rip.pyr.Rall=[sess.rip.pyr.Rall;xml];
      end
    end
    %% Pyr Only, Ripple Mod vs Nomod
    if enough.rippledata
      if enough.pyr.Rmod % Minimum 6 cells in each structure
	EVREV.sws.pyr.Rmod=[EVREV.sws.pyr.Rmod; EV.sws.pyr.Rmod REV.sws.pyr.Rmod EV.post1.sws.pyr.Rmod REV.post1.sws.pyr.Rmod EV.post2.sws.pyr.Rmod REV.post2.sws.pyr.Rmod EV.post3.sws.pyr.Rmod REV.post3.sws.pyr.Rmod];
	EVREV.sws.pyr.Rnomod=[EVREV.sws.pyr.Rnomod; EV.sws.pyr.Rnomod REV.sws.pyr.Rnomod EV.post1.sws.pyr.Rnomod REV.post1.sws.pyr.Rnomod EV.post2.sws.pyr.Rnomod REV.post2.sws.pyr.Rnomod EV.post3.sws.pyr.Rnomod REV.post3.sws.pyr.Rnomod];
	sess.sws.pyr.Rmod=[sess.sws.pyr.Rmod;xml];
	if enough.rem
	  EVREV.rem.pyr.Rmod=[EVREV.rem.pyr.Rmod; EV.rem.pyr.Rmod REV.rem.pyr.Rmod];
	  EVREV.rem.pyr.Rnomod=[EVREV.rem.pyr.Rnomod; EV.rem.pyr.Rnomod REV.rem.pyr.Rnomod];
	  sess.rem.pyr.Rmod=[sess.rem.pyr.Rmod;xml];
	end
	if treatrun
	  EVREV.run.pyr.Rmod=[EVREV.run.pyr.Rmod; EV.run.pyr.Rmod REV.run.pyr.Rmod];
	  EVREV.run.pyr.Rnomod=[EVREV.run.pyr.Rnomod; EV.run.pyr.Rnomod REV.run.pyr.Rnomod];
	  sess.run.pyr.Rmod=[sess.run.pyr.Rmod;xml];
	end
	EVREV.rip_in.pyr.Rmod=[EVREV.rip_in.pyr.Rmod; EV.rip_in.pyr.Rmod REV.rip_in.pyr.Rmod];
	EVREV.rip_in.pyr.Rnomod=[EVREV.rip_in.pyr.Rnomod; EV.rip_in.pyr.Rnomod REV.rip_in.pyr.Rnomod];
	EVREV.rip_out.pyr.Rmod=[EVREV.rip_out.pyr.Rmod; EV.rip_out.pyr.Rmod REV.rip_out.pyr.Rmod];
	EVREV.rip_out.pyr.Rnomod=[EVREV.rip_out.pyr.Rnomod; EV.rip_out.pyr.Rnomod REV.rip_out.pyr.Rnomod];
	sess.rip.pyr.Rmod=[sess.rip.pyr.Rmod;xml];
      end
    end
  end
end
    
%%%%%%%%%

meanEVREV.sws.all.Rall=nanmean(EVREV.sws.all.Rall);
swsEVREV.all.Rall=EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),:);
swsmeanEVREV.all.Rall=nanmean(swsEVREV.all.Rall);
meanEVREV.rem.all.Rall=mean(EVREV.rem.all.Rall);
meanEVREV.rip_in.all.Rall=mean(EVREV.rip_in.all.Rall);
meanEVREV.rip_out.all.Rall=mean(EVREV.rip_out.all.Rall);

meanEVREV.sws.pyr.Rall=nanmean(EVREV.sws.pyr.Rall);
swsEVREV.pyr.Rall=EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),:);
swsmeanEVREV.pyr.Rall=nanmean(swsEVREV.pyr.Rall);
meanEVREV.rem.pyr.Rall=mean(EVREV.rem.pyr.Rall);
meanEVREV.rip_in.pyr.Rall=mean(EVREV.rip_in.pyr.Rall);
meanEVREV.rip_out.pyr.Rall=mean(EVREV.rip_out.pyr.Rall);

meanEVREV.run.all.Rall=nanmean(EVREV.run.all.Rall);
meanEVREV.run.pyr.Rall=nanmean(EVREV.run.pyr.Rall);

semEVREV.sws.all.Rall=nansem(EVREV.sws.all.Rall);
swssemEVREV.all.Rall=nansem(swsEVREV.all.Rall);
semEVREV.rem.all.Rall=sem(EVREV.rem.all.Rall);
semEVREV.rip_in.all.Rall=sem(EVREV.rip_in.all.Rall);
semEVREV.rip_out.all.Rall=sem(EVREV.rip_out.all.Rall);
semEVREV.sws.pyr.Rall=nansem(EVREV.sws.pyr.Rall);
swssemEVREV.pyr.Rall=nansem(swsEVREV.pyr.Rall);
semEVREV.rem.pyr.Rall=sem(EVREV.rem.pyr.Rall);
semEVREV.rip_in.pyr.Rall=sem(EVREV.rip_in.pyr.Rall);
semEVREV.rip_out.pyr.Rall=sem(EVREV.rip_out.pyr.Rall);

%  keyboard
if treatrun
  semEVREV.run.all.Rall=nansem(EVREV.run.all.Rall);
  semEVREV.run.pyr.Rall=nansem(EVREV.run.pyr.Rall);
end

intrahpc=strcmp(structure,'Hpc') & strcmp(EVtype,'intra');

if skiprip
  disp('No data for ripple modualtion - intra Hpc')
else
  meanEVREV.rem.all.Rmod=mean(EVREV.rem.all.Rmod);
  meanEVREV.rem.all.Rnomod=mean(EVREV.rem.all.Rnomod);
  meanEVREV.rip_in.all.Rmod=mean(EVREV.rip_in.all.Rmod);
  meanEVREV.rip_in.all.Rnomod=mean(EVREV.rip_in.all.Rnomod);
  meanEVREV.rip_out.all.Rmod=mean(EVREV.rip_out.all.Rmod);
  meanEVREV.rip_out.all.Rnomod=mean(EVREV.rip_out.all.Rnomod);
  meanEVREV.sws.all.Rmod=nanmean(EVREV.sws.all.Rmod);
  meanEVREV.sws.all.Rnomod=nanmean(EVREV.sws.all.Rnomod);
  swsmeanEVREV.all.Rmod=nanmean(EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),:));
  swsmeanEVREV.all.Rnomod=nanmean(EVREV.sws.all.Rnomod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),:));
  
%    meanEVREV.run.all.Rmod=nanmean(EVREV.run.all.Rmod);
%    meanEVREV.run.all.Rnomod=nanmean(EVREV.run.all.Rnomod);
%    semEVREV.run.all.Rmod=nansem(EVREV.run.all.Rmod);
%    semEVREV.run.all.Rnomod=nansem(EVREV.run.all.Rnomod);
%    
%    meanEVREV.run.pyr.Rmod=nanmean(EVREV.run.pyr.Rmod);
%    meanEVREV.run.pyr.Rnomod=nanmean(EVREV.run.pyr.Rnomod);
%    semEVREV.run.pyr.Rmod=nansem(EVREV.run.pyr.Rmod);
%    semEVREV.run.pyr.Rnomod=nansem(EVREV.run.pyr.Rnomod);
  
  meanEVREV.sws.pyr.Rmod=nanmean(EVREV.sws.pyr.Rmod);
  meanEVREV.sws.pyr.Rnomod=nanmean(EVREV.sws.pyr.Rnomod);
  swsmeanEVREV.pyr.Rmod=nanmean(EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),:));
  swsmeanEVREV.pyr.Rnomod=nanmean(EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),:));

  meanEVREV.rem.pyr.Rmod=mean(EVREV.rem.pyr.Rmod);
  meanEVREV.rem.pyr.Rnomod=mean(EVREV.rem.pyr.Rnomod);
  meanEVREV.rip_in.pyr.Rmod=mean(EVREV.rip_in.pyr.Rmod);
  meanEVREV.rip_in.pyr.Rnomod=mean(EVREV.rip_in.pyr.Rnomod);
  meanEVREV.rip_out.pyr.Rmod=mean(EVREV.rip_out.pyr.Rmod);
  meanEVREV.rip_out.pyr.Rnomod=mean(EVREV.rip_out.pyr.Rnomod);

  semEVREV.sws.all.Rmod=nansem(EVREV.sws.all.Rmod);
  semEVREV.sws.all.Rnomod=nansem(EVREV.sws.all.Rnomod);
  swssemEVREV.all.Rmod=nansem(EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),:));
  swssemEVREV.all.Rnomod=nansem(EVREV.sws.all.Rnomod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),:));
  semEVREV.rem.all.Rmod=sem(EVREV.rem.all.Rmod);
  semEVREV.rem.all.Rnomod=sem(EVREV.rem.all.Rnomod);
  semEVREV.rip_in.all.Rmod=sem(EVREV.rip_in.all.Rmod);
  semEVREV.rip_in.all.Rnomod=sem(EVREV.rip_in.all.Rnomod);
  semEVREV.rip_out.all.Rmod=sem(EVREV.rip_out.all.Rmod);
  semEVREV.rip_out.all.Rnomod=sem(EVREV.rip_out.all.Rnomod);
  semEVREV.sws.pyr.Rmod=nansem(EVREV.sws.pyr.Rmod);
  semEVREV.sws.pyr.Rnomod=nansem(EVREV.sws.pyr.Rnomod);
  
  swssemEVREV.pyr.Rmod=nansem(EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),:));
  swssemEVREV.pyr.Rnomod=nansem(EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),:));
  semEVREV.rem.pyr.Rmod=sem(EVREV.rem.pyr.Rmod);
  semEVREV.rem.pyr.Rnomod=sem(EVREV.rem.pyr.Rnomod);
  semEVREV.rip_in.pyr.Rmod=sem(EVREV.rip_in.pyr.Rmod);
  semEVREV.rip_in.pyr.Rnomod=sem(EVREV.rip_in.pyr.Rnomod);
  semEVREV.rip_out.pyr.Rmod=sem(EVREV.rip_out.pyr.Rmod);
  semEVREV.rip_out.pyr.Rnomod=sem(EVREV.rip_out.pyr.Rnomod);
end


%%%%% Figures
titsws={'EV','REV','','EV-1','REV-1','EV-2','REV2','EV-3','REV-3'};
titsws1={'EV','REV','','EV-1','REV-1','EV-2','REV2','EV-3','REV-3'};
titrem={'EV','REV'};
titrun={'EV','REV'};
%  titswsMod={'EV','','REV','','','EV-1','','REV-1','','EV-2','','REV2','','EV-3','','REV-3'};
if strcmp(EVtype,'cross')
  suptit=['HPC-' struct ' pairs'];
elseif strcmp(EVtype,'intra')
  suptit=['Intra-' struct ' pairs'];
end

%  %Stats figure 1
[p.sws.all.Rall,~,stats.sws.all.Rall]=signrank(EVREV.sws.all.Rall(:,1),EVREV.sws.all.Rall(:,2));
[p.post1.all.Rall,~,stats.post1.all.Rall]=signrank(EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),3),EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),4));
[p.post2.all.Rall,~,stats.post2.all.Rall]=signrank(EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),5),EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),6));
[p.post3.all.Rall,~,stats.post3.all.Rall]=signrank(EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),7),EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),8));

[p.sws.pyr.Rall,~,stats.sws.pyr.Rall]=signrank(EVREV.sws.pyr.Rall(:,1),EVREV.sws.pyr.Rall(:,2));
[p.post1.pyr.Rall,~,stats.post1.pyr.Rall]=signrank(EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),3),EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),4));
[p.post2.pyr.Rall,~,stats.post2.pyr.Rall]=signrank(EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),5),EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),6));
[p.post3.pyr.Rall,~,stats.post3.pyr.Rall]=signrank(EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),7),EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),8));

if ~skiprip
  [p.sws.all.Rmod,~,stats.sws.all.Rmod]=signrank(EVREV.sws.all.Rmod(:,1),EVREV.sws.all.Rmod(:,2));
  p.post1.all.Rmod=signrank(EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),3),EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),4));
  p.post2.all.Rmod=signrank(EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),5),EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),6));
  p.post3.all.Rmod=signrank(EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),7),EVREV.sws.all.Rmod(EVREV.sws.all.Rmod(:,3)>EVREV.sws.all.Rmod(:,4),8));

  [p.sws.all.Rnomod,~,stats.sws.all.Rnomod]=signrank(EVREV.sws.all.Rnomod(:,1),EVREV.sws.all.Rnomod(:,2));
  p.post1.all.Rnomod=signrank(EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),3),EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),4));
  p.post2.all.Rnomod=signrank(EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),5),EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),6));
  p.post3.all.Rnomod=signrank(EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),7),EVREV.sws.all.Rnomod(EVREV.sws.all.Rnomod(:,3)>EVREV.sws.all.Rnomod(:,4),8));

  [p.sws.pyr.Rmod,~,stats.sws.pyr.Rmod]=signrank(EVREV.sws.pyr.Rmod(:,1),EVREV.sws.pyr.Rmod(:,2));
  p.post1.pyr.Rmod=signrank(EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),3),EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),4));
  p.post2.pyr.Rmod=signrank(EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),5),EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),6));
  p.post3.pyr.Rmod=signrank(EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),7),EVREV.sws.pyr.Rmod(EVREV.sws.pyr.Rmod(:,3)>EVREV.sws.pyr.Rmod(:,4),8));

  [p.sws.pyr.Rnomod,~,stats.sws.pyr.Rnomod]=signrank(EVREV.sws.pyr.Rnomod(:,1),EVREV.sws.pyr.Rnomod(:,2));
  p.post1.pyr.Rnomod=signrank(EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),3),EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),4));
  p.post2.pyr.Rnomod=signrank(EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),5),EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),6));
  p.post3.pyr.Rnomod=signrank(EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),7),EVREV.sws.pyr.Rnomod(EVREV.sws.pyr.Rnomod(:,3)>EVREV.sws.pyr.Rnomod(:,4),8));
end

%% EVREV differences and test
diff.all.post1=EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),3)-EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),4);
diff.all.post2=EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),5)-EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),6);
diff.all.post3=EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),7)-EVREV.sws.all.Rall(EVREV.sws.all.Rall(:,3)>EVREV.sws.all.Rall(:,4),8);

diff.pyr.post1=EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),3)-EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),4);
diff.pyr.post2=EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),5)-EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),6);
diff.pyr.post3=EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),7)-EVREV.sws.pyr.Rall(EVREV.sws.pyr.Rall(:,3)>EVREV.sws.pyr.Rall(:,4),8);

diff.all.Rmod=EVREV.sws.all.Rmod(:,1)-EVREV.sws.all.Rmod(:,2);
diff.all.Rnomod=EVREV.sws.all.Rnomod(:,1)-EVREV.sws.all.Rnomod(:,2);
diff.pyr.Rmod=EVREV.sws.pyr.Rmod(:,1)-EVREV.sws.pyr.Rmod(:,2);
diff.pyr.Rnomod=EVREV.sws.pyr.Rnomod(:,1)-EVREV.sws.pyr.Rnomod(:,2);

%  % Test for diff ripple mod/nomod
%  [p,h,st]=ranksum(diff.all.Rmod,diff.all.Rnomod)
%  [p,h,st]=ranksum(diff.pyr.Rmod,diff.pyr.Rnomod)
%  
%  
%  % Test for decay, pyr pairs
%  [pkw,tab,sts]=kruskalwallis([diff.pyr.post1 diff.pyr.post2 diff.pyr.post3]);
%  [c,m,h,nms]=multcompare(sts);

% BARPLOT FIGURES
%  %%
%  f1=figure('Position',[370 76 1280 877]);
%  subplot(2,2,1);
%  barwitherr([semEVREV.sws.all.Rall(:,1:2) 0 swssemEVREV.all.Rall(:,3:8)],[meanEVREV.sws.all.Rall(:,1:2) 0  swsmeanEVREV.all.Rall(:,3:8)],'FaceColor',rgb('Crimson'),'EdgeColor','none');
%  set(gca,'XTickLabel',titsws);
%  xlabel('SWS - all cell types');
%  ylabel(['n=' int2str(length(sess.sws.all.Rall)) 'sess']);
%  sigstar({[1 2],[4 5],[6 7],[8 9]},[p.sws.all.Rall,p.post1.all.Rall,p.post2.all.Rall,p.post3.all.Rall]);
%  hold on;
%  plot([EVREV.sws.all.Rall(:,1:2)]','-k');
%  plot([4 5],swsEVREV.all.Rall(:,3:4)','-k');
%  plot([6 7],swsEVREV.all.Rall(:,5:6)','-k');
%  plot([8 9],swsEVREV.all.Rall(:,7:8)','-k');
%  %  %  keyboard
%  
%  if ~skiprip
%    subplot(2,2,2);
%    barwitherr([semEVREV.sws.all.Rmod(:,1) semEVREV.sws.all.Rmod(:,2) 0 0 0 swssemEVREV.all.Rmod(:,3) 0 swssemEVREV.all.Rmod(:,4) 0 swssemEVREV.all.Rmod(:,5) 0 swssemEVREV.all.Rmod(:,6) 0 swssemEVREV.all.Rmod(:,7) 0 swssemEVREV.all.Rmod(:,8) 0],[meanEVREV.sws.all.Rmod(:,1) meanEVREV.sws.all.Rmod(:,2) 0 0 0 swsmeanEVREV.all.Rmod(:,3) 0 swsmeanEVREV.all.Rmod(:,4) 0 swsmeanEVREV.all.Rmod(:,5) 0 swsmeanEVREV.all.Rmod(:,6) 0 swsmeanEVREV.all.Rmod(:,7) 0 swsmeanEVREV.all.Rmod(:,8) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
%    %  set(gca,'XTickLabel',titsws1);
%    hold on;
%    barwitherr([0 0 semEVREV.sws.all.Rnomod(:,1) semEVREV.sws.all.Rnomod(:,2) 0 0 swssemEVREV.all.Rnomod(:,3) 0 swssemEVREV.all.Rnomod(:,4) 0 swssemEVREV.all.Rnomod(:,5) 0 swssemEVREV.all.Rnomod(:,6) 0 swssemEVREV.all.Rnomod(:,7) 0 swssemEVREV.all.Rnomod(:,8)],[0 0 meanEVREV.sws.all.Rnomod(:,1) meanEVREV.sws.all.Rnomod(:,2) 0 0 swsmeanEVREV.all.Rnomod(:,3) 0 swsmeanEVREV.all.Rnomod(:,4) 0 swsmeanEVREV.all.Rnomod(:,5) 0 swsmeanEVREV.all.Rnomod(:,6) 0 swsmeanEVREV.all.Rnomod(:,7) 0 swsmeanEVREV.all.Rnomod(:,8)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
%    xlabel('Ripple-modulated cells (Red) vs Non-rip Mod cells (Green) - All cell types');
%    ylabel(['n=' int2str(length(sess.sws.all.Rmod)) 'sess']);
%    sigstar({[1 3],[2 4],[6 8],[7 9],[10 12],[11 13],[14 16],[15 17]},[p.sws.all.Rmod,p.sws.all.Rnomod,p.post1.all.Rmod,p.post1.all.Rnomod,p.post2.all.Rmod,p.post2.all.Rnomod,p.post3.all.Rmod,p.post3.all.Rnomod]);
%    hold on;
%    plot(EVREV.sws.all.Rmod(:,1:2)','-k');
%    plot([3 4],EVREV.sws.all.Rnomod(:,1:2)','-k');
%    
%    subplot(2,2,4);
%    barwitherr([semEVREV.sws.pyr.Rmod(:,1) semEVREV.sws.pyr.Rmod(:,2) 0 0 0 swssemEVREV.pyr.Rmod(:,3) 0 swssemEVREV.pyr.Rmod(:,4) 0 swssemEVREV.pyr.Rmod(:,5) 0 swssemEVREV.pyr.Rmod(:,6) 0 swssemEVREV.pyr.Rmod(:,7) 0 swssemEVREV.pyr.Rmod(:,8) 0],[meanEVREV.sws.pyr.Rmod(:,1) meanEVREV.sws.pyr.Rmod(:,2) 0 0 0 swsmeanEVREV.pyr.Rmod(:,3) 0 swsmeanEVREV.pyr.Rmod(:,4) 0 swsmeanEVREV.pyr.Rmod(:,5) 0 swsmeanEVREV.pyr.Rmod(:,6) 0 swsmeanEVREV.pyr.Rmod(:,7) 0 swsmeanEVREV.pyr.Rmod(:,8) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
%    %  set(gca,'XTickLabel',titsws1);
%    hold on;
%    barwitherr([0 0 semEVREV.sws.pyr.Rnomod(:,1) semEVREV.sws.pyr.Rnomod(:,2) 0 0 swssemEVREV.pyr.Rnomod(:,3) 0 swssemEVREV.pyr.Rnomod(:,4) 0 swssemEVREV.pyr.Rnomod(:,5) 0 swssemEVREV.pyr.Rnomod(:,6) 0 swssemEVREV.pyr.Rnomod(:,7) 0 swssemEVREV.pyr.Rnomod(:,8)],[0 0 meanEVREV.sws.pyr.Rnomod(:,1) meanEVREV.sws.pyr.Rnomod(:,2) 0 0 swsmeanEVREV.pyr.Rnomod(:,3) 0 swsmeanEVREV.pyr.Rnomod(:,4) 0 swsmeanEVREV.pyr.Rnomod(:,5) 0 swsmeanEVREV.pyr.Rnomod(:,6) 0 swsmeanEVREV.pyr.Rnomod(:,7) 0 swsmeanEVREV.pyr.Rnomod(:,8)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
%    xlabel('Ripple-modulated cells (REd) vs Non-rip Mod cells (Green) - Pyr Only');
%    ylabel(['n=' int2str(length(sess.sws.pyr.Rmod)) 'sess']);
%    sigstar({[1 3],[2 4],[6 8],[7 9],[10 12],[11 13],[14 16],[15 17]},[p.sws.pyr.Rmod,p.sws.pyr.Rnomod,p.post1.pyr.Rmod,p.post1.pyr.Rnomod,p.post2.pyr.Rmod,p.post2.pyr.Rnomod,p.post3.pyr.Rmod,p.post3.pyr.Rnomod]);
%    hold on;
%    plot(EVREV.sws.pyr.Rmod(:,1:2)','-k');
%    plot([3 4],EVREV.sws.pyr.Rnomod(:,1:2)','-k');  
%  end
%  
%  subplot(2,2,3);
%  barwitherr([semEVREV.sws.pyr.Rall(:,1:2) 0 swssemEVREV.pyr.Rall(:,3:8)],[meanEVREV.sws.pyr.Rall(:,1:2) 0  swsmeanEVREV.pyr.Rall(:,3:8)],'FaceColor',rgb('Crimson'),'EdgeColor','none');
%  set(gca,'XTickLabel',titsws);
%  xlabel('SWS - Pyr Only');
%  ylabel(['n=' int2str(length(sess.sws.pyr.Rall)) 'sess']);
%  sigstar({[1 2],[4 5],[6 7],[8 9]},[p.sws.pyr.Rall,p.post1.pyr.Rall,p.post2.pyr.Rall,p.post3.pyr.Rall]);
%  hold on;
%  plot([EVREV.sws.pyr.Rall(:,1:2)]','-k');
%  plot([4 5],swsEVREV.pyr.Rall(:,3:4)','-k');
%  plot([6 7],swsEVREV.pyr.Rall(:,5:6)','-k');
%  plot([8 9],swsEVREV.pyr.Rall(:,7:8)','-k');
%  suptitle(suptit);

%%%%%%%%%%%%%%%%% Boxplot Figs    
f1=figure('Position',[370 76 1280 877]);
subplot(2,2,1);
BPswsEVREV.all.Rall=swsEVREV.all.Rall(:,3:8);
d1=[EVREV.sws.all.Rall(:,1);EVREV.sws.all.Rall(:,2);BPswsEVREV.all.Rall(:)]
a=repmat([1 2],size(EVREV.sws.all.Rall,1),1);a=a(:);
b=repmat([3:8],size(BPswsEVREV.all.Rall,1),1);b=b(:);
g1=[a;b];
c1=mod(g1,2);
boxplot(d1,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);
set(gca,'XTickLabel',titsws);
xlabel('SWS - all cell types');
ylabel(['n=' int2str(length(sess.sws.all.Rall)) 'sess']);
sigstar({[1 2],[3 4],[5 6],[7 8]},[p.sws.all.Rall,p.post1.all.Rall,p.post2.all.Rall,p.post3.all.Rall]);

subplot(2,2,3);
BPswsEVREV.pyr.Rall=swsEVREV.pyr.Rall(:,3:8);
d1=[EVREV.sws.pyr.Rall(:,1);EVREV.sws.pyr.Rall(:,2);BPswsEVREV.pyr.Rall(:)]
a=repmat([1 2],size(EVREV.sws.pyr.Rall,1),1);a=a(:);
b=repmat([3:8],size(BPswsEVREV.pyr.Rall,1),1);b=b(:);
g1=[a;b];
c1=mod(g1,2);
boxplot(d1,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);
set(gca,'XTickLabel',titsws);
xlabel('SWS - Pyr Only');
ylabel(['n=' int2str(length(sess.sws.pyr.Rall)) 'sess']);
sigstar({[1 2],[3 4],[5 6],[7 8]},[p.sws.pyr.Rall,p.post1.pyr.Rall,p.post2.pyr.Rall,p.post3.pyr.Rall]);

if ~skiprip
  subplot(2,2,2);
  d1=[EVREV.sws.all.Rmod(:,1);EVREV.sws.all.Rmod(:,2);EVREV.sws.all.Rnomod(:,1);EVREV.sws.all.Rnomod(:,2)];
  a=repmat([1 2],size(EVREV.sws.all.Rmod,1),1);a=a(:);
  b=repmat([3 4],size(EVREV.sws.all.Rnomod,1),1);b=b(:);
  g1=[a;b];
  c1=mod(g1,2);%odd/even for colors
  boxplot(d1,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);
  ylabel(['n=' int2str(length(sess.sws.all.Rmod)) 'sess']);
  xlabel('Ripple-modulated cells (REd) vs Non-rip Mod cells (Green) - All Cells');
  sigstar({[1 2],[3 4]},[p.sws.all.Rmod,p.sws.all.Rnomod]);
  
  subplot(2,2,4);
  d1=[EVREV.sws.pyr.Rmod(:,1);EVREV.sws.pyr.Rmod(:,2);EVREV.sws.pyr.Rnomod(:,1);EVREV.sws.pyr.Rnomod(:,2)];
  a=repmat([1 2],size(EVREV.sws.pyr.Rmod,1),1);a=a(:);
  b=repmat([3 4],size(EVREV.sws.pyr.Rnomod,1),1);b=b(:);
  g1=[a;b];
  c1=mod(g1,2);%odd/even for colors
  boxplot(d1,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);
  xlabel('Ripple-modulated cells (REd) vs Non-rip Mod cells (Green) - Pyr Only');
  ylabel(['n=' int2str(length(sess.sws.pyr.Rmod)) 'sess']);
  sigstar({[1 2],[3 4]},[p.sws.pyr.Rmod,p.sws.pyr.Rnomod]);  
end


% Stats figure2
[p.rem.all.Rall,~,stats.rem.all.Rall]=signrank(EVREV.rem.all.Rall(:,1),EVREV.rem.all.Rall(:,2));
[p.rem.pyr.Rall,~,stats.rem.pyr.Rall]=signrank(EVREV.rem.pyr.Rall(:,1),EVREV.rem.pyr.Rall(:,2));
if ~skiprip
  p.rem.all.Rmod=signrank(EVREV.rem.all.Rmod(:,1),EVREV.rem.all.Rmod(:,2));
  p.rem.all.Rnomod=signrank(EVREV.rem.all.Rnomod(:,1),EVREV.rem.all.Rnomod(:,2));
  p.rem.pyr.Rmod=signrank(EVREV.rem.pyr.Rmod(:,1),EVREV.rem.pyr.Rmod(:,2));
  p.rem.pyr.Rnomod=signrank(EVREV.rem.pyr.Rnomod(:,1),EVREV.rem.pyr.Rnomod(:,2));
end


%%
f2=figure('Position',[1337 111 401 835]);
subplot(2,2,1);
%  barwitherr(semEVREV.rem.all.Rall(:,1:2),meanEVREV.rem.all.Rall(:,1:2),'FaceColor',rgb('Crimson'),'EdgeColor','none');
%  hold on;
%  plot([1 2],EVREV.rem.all.Rall(:,1:2)','-k');
g1=repmat([1 2],size(EVREV.rem.all.Rall,1),1);g1=g1(:);
c1=mod(g1,2);
boxplot(EVREV.rem.all.Rall,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);
set(gca,'XTickLabel',titrem);
sigstar([1 2],p.rem.all.Rall);
xlabel('REM All cells');
ylabel(['n=' int2str(length(sess.rem.all.Rall)) 'sess']);

if ~skiprip
  subplot(2,2,2);
  barwitherr([semEVREV.rem.all.Rmod(:,1) 0 semEVREV.rem.all.Rmod(:,2) 0],[meanEVREV.rem.all.Rmod(:,1) 0 meanEVREV.rem.all.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
  hold on;
  barwitherr([0 semEVREV.rem.all.Rnomod(:,1) 0 semEVREV.rem.all.Rnomod(:,2)],[0 meanEVREV.rem.all.Rnomod(:,1) 0 meanEVREV.rem.all.Rnomod(:,2)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
  xlabel('REM All cells Rip Mod/NoMod');
  ylabel(['n=' int2str(length(sess.rem.all.Rmod)) 'sess']);
  sigstar({[1 3],[2 4]},[p.rem.all.Rmod,p.rem.all.Rnomod]);
  subplot(2,2,4);
  barwitherr([semEVREV.rem.pyr.Rmod(:,1) 0 semEVREV.rem.pyr.Rmod(:,2) 0],[meanEVREV.rem.pyr.Rmod(:,1) 0 meanEVREV.rem.pyr.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
  hold on;
  barwitherr([0 semEVREV.rem.pyr.Rnomod(:,1) 0 semEVREV.rem.pyr.Rnomod(:,2)],[0 meanEVREV.rem.pyr.Rnomod(:,1) 0 meanEVREV.rem.pyr.Rnomod(:,2)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
  xlabel('REM Pyr Only Rip Mod/NoMod');
  ylabel(['n=' int2str(length(sess.rem.pyr.Rmod)) 'sess']);
  sigstar({[1 3],[2 4]},[p.rem.pyr.Rmod,p.rem.pyr.Rnomod]);
end
subplot(2,2,3);
%  barwitherr(semEVREV.rem.pyr.Rall(:,1:2),meanEVREV.rem.pyr.Rall(:,1:2),'FaceColor',rgb('Crimson'),'EdgeColor','none');
%  hold on;
%  plot([1 2],EVREV.rem.pyr.Rall(:,1:2)','-k');
g1=repmat([1 2],size(EVREV.rem.pyr.Rall,1),1);g1=g1(:);
c1=mod(g1,2);
boxplot(EVREV.rem.pyr.Rall,g1,'MedianStyle','line','Symbol','k.','OutlierSize',8,'Jitter',0,'Widths',0.2,'ColorGroup',c1);

set(gca,'XTickLabel',titrem);
xlabel('REM Pyr Only');
ylabel(['n=' int2str(length(sess.rem.pyr.Rall)) 'sess']);
sigstar([1 2],p.rem.pyr.Rall);
suptitle(suptit);


%  % Stats Figure3
p.rip_in.all.Rall=signrank(EVREV.rip_in.all.Rall(:,1),EVREV.rip_in.all.Rall(:,2));
p.rip_in.pyr.Rall=signrank(EVREV.rip_in.pyr.Rall(:,1),EVREV.rip_in.pyr.Rall(:,2));
p.rip_out.all.Rall=signrank(EVREV.rip_out.all.Rall(:,1),EVREV.rip_out.all.Rall(:,2));
p.rip_out.pyr.Rall=signrank(EVREV.rip_out.pyr.Rall(:,1),EVREV.rip_out.pyr.Rall(:,2));

if ~skiprip
  p.rip_in.all.Rmod=signrank(EVREV.rip_in.all.Rmod(:,1),EVREV.rip_in.all.Rmod(:,2));
  p.rip_in.all.Rnomod=signrank(EVREV.rip_in.all.Rnomod(:,1),EVREV.rip_in.all.Rnomod(:,2));
  p.rip_in.pyr.Rmod=signrank(EVREV.rip_in.pyr.Rmod(:,1),EVREV.rip_in.pyr.Rmod(:,2));
  p.rip_in.pyr.Rnomod=signrank(EVREV.rip_in.pyr.Rnomod(:,1),EVREV.rip_in.pyr.Rnomod(:,2));
  p.rip_out.all.Rmod=signrank(EVREV.rip_out.all.Rmod(:,1),EVREV.rip_out.all.Rmod(:,2));
  p.rip_out.all.Rnomod=signrank(EVREV.rip_out.all.Rnomod(:,1),EVREV.rip_out.all.Rnomod(:,2));
  p.rip_out.pyr.Rmod=signrank(EVREV.rip_out.pyr.Rmod(:,1),EVREV.rip_out.pyr.Rmod(:,2));
  p.rip_out.pyr.Rnomod=signrank(EVREV.rip_out.pyr.Rnomod(:,1),EVREV.rip_out.pyr.Rnomod(:,2));
end

%%
f3=figure('Position',[781 113 549 835]);
subplot(2,2,1);
barwitherr([semEVREV.rip_in.all.Rall(:,1) 0 semEVREV.rip_in.all.Rall(:,2) 0],[meanEVREV.rip_in.all.Rall(:,1) 0 meanEVREV.rip_in.all.Rall(2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
hold on;
barwitherr([0 semEVREV.rip_out.all.Rall(:,1) 0 semEVREV.rip_out.all.Rall(:,2)],[0 meanEVREV.rip_out.all.Rall(:,1) 0 meanEVREV.rip_out.all.Rall(2)],'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
sigstar({[1 3],[2 4]},[p.rip_in.all.Rall,p.rip_out.all.Rall]);
xlabel('All cells');
ylabel(['n=' int2str(length(sess.rip.all.Rall)) 'sess']);
set(gca,'XTickLabel',{'EV' 'EV' 'REV' 'REV'});
if ~skiprip
  subplot(2,2,2);
  barwitherr([semEVREV.rip_in.all.Rmod(:,1) 0 semEVREV.rip_in.all.Rmod(:,2) 0],[meanEVREV.rip_in.all.Rmod(:,1) 0 meanEVREV.rip_in.all.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
  hold on;
  barwitherr([0 semEVREV.rip_out.all.Rmod(:,1) 0 semEVREV.rip_out.all.Rmod(:,2)],[0 meanEVREV.rip_out.all.Rmod(:,1) 0 meanEVREV.rip_out.all.Rmod(:,2)],'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
  xlabel('All Ripple Mod cells');
  ylabel(['n=' int2str(length(sess.rip.all.Rmod)) 'sess']);
  set(gca,'XTickLabel',{'EV' 'EV' 'REV' 'REV'});
  sigstar({[1 3],[2 4]},[p.rip_in.all.Rmod,p.rip_out.all.Rmod]);
  subplot(2,2,4);
  barwitherr([semEVREV.rip_in.pyr.Rmod(:,1) 0 semEVREV.rip_in.pyr.Rmod(:,2) 0],[meanEVREV.rip_in.pyr.Rmod(:,1) 0 meanEVREV.rip_in.pyr.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
  hold on;
  barwitherr([0 semEVREV.rip_out.pyr.Rmod(:,1) 0 semEVREV.rip_out.pyr.Rmod(:,2)],[0 meanEVREV.rip_out.pyr.Rmod(:,1) 0 meanEVREV.rip_out.pyr.Rmod(:,2)],'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
  xlabel('Mod Pyr');
  ylabel(['n=' int2str(length(sess.rip.pyr.Rmod)) 'sess']);
  sigstar({[1 3],[2 4]},[p.rip_in.pyr.Rmod,p.rip_out.pyr.Rmod]);
  set(gca,'XTickLabel',{'EV' 'EV' 'REV' 'REV'});
end
subplot(2,2,3);
barwitherr([semEVREV.rip_in.pyr.Rall(1) 0 semEVREV.rip_in.pyr.Rall(2) 0],[meanEVREV.rip_in.pyr.Rall(1) 0 meanEVREV.rip_in.pyr.Rall(2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
hold on
barwitherr([0 semEVREV.rip_out.pyr.Rall(1) 0 semEVREV.rip_in.pyr.Rall(2)],[0 meanEVREV.rip_out.pyr.Rall(1) 0 meanEVREV.rip_out.pyr.Rall(2)],'FaceColor',rgb('MidnightBlue'),'EdgeColor','none');
xlabel('Pyr Only');
ylabel(['n=' int2str(length(sess.rip.pyr.Rall)) 'sess']);
sigstar({[1 3],[2 4]},[p.rip_in.pyr.Rall,p.rip_out.pyr.Rall]);
set(gca,'XTickLabel',{'EV' 'EV' 'REV' 'REV'});
suptitle([suptit '-In (Red) vs Out (Blue) Ripples']);    

%% Stats figure 4
if treatrun
  p.run.all.Rall=signrank(EVREV.run.all.Rall(:,1),EVREV.run.all.Rall(:,2));
  p.run.pyr.Rall=signrank(EVREV.run.pyr.Rall(:,1),EVREV.run.pyr.Rall(:,2));
  if ~skiprip
    p.run.all.Rmod=signrank(EVREV.run.all.Rmod(:,1),EVREV.run.all.Rmod(:,2));
    p.run.all.Rnomod=signrank(EVREV.run.all.Rnomod(:,1),EVREV.run.all.Rnomod(:,2));
    p.run.pyr.Rmod=signrank(EVREV.run.pyr.Rmod(:,1),EVREV.run.pyr.Rmod(:,2));
    p.run.pyr.Rnomod=signrank(EVREV.run.pyr.Rnomod(:,1),EVREV.run.pyr.Rnomod(:,2));
  end
end

if treatrun
  f4=figure('Position',[1337 111 401 835]);
  subplot(2,2,1);
  barwitherr(semEVREV.run.all.Rall(:,1:2),meanEVREV.run.all.Rall(:,1:2),'FaceColor',rgb('Crimson'),'EdgeColor','none');
  set(gca,'XTickLabel',titrun);
  xlabel('RUN All cells');
  ylabel(['n=' int2str(length(sess.run.all.Rall)) 'sess']);
  sigstar([1 2],p.run.all.Rall);
  if ~skiprip
    subplot(2,2,2);
    barwitherr([semEVREV.run.all.Rmod(:,1) 0 semEVREV.run.all.Rmod(:,2) 0],[meanEVREV.run.all.Rmod(:,1) 0 meanEVREV.run.all.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
    hold on;
    barwitherr([0 semEVREV.run.all.Rnomod(:,1) 0 semEVREV.run.all.Rnomod(:,2)],[0 meanEVREV.run.all.Rnomod(:,1) 0 meanEVREV.run.all.Rnomod(:,2)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
    xlabel('RUN All cells Rip Mod/NoMod');
    ylabel(['n=' int2str(length(sess.run.all.Rmod)) 'sess']);
    sigstar({[1 3],[2 4]},[p.run.all.Rmod,p.run.all.Rnomod]);
    subplot(2,2,4);
    barwitherr([semEVREV.run.pyr.Rmod(:,1) 0 semEVREV.run.pyr.Rmod(:,2) 0],[meanEVREV.run.pyr.Rmod(:,1) 0 meanEVREV.run.pyr.Rmod(:,2) 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');
    hold on;
    barwitherr([0 semEVREV.run.pyr.Rnomod(:,1) 0 semEVREV.run.pyr.Rnomod(:,2)],[0 meanEVREV.run.pyr.Rnomod(:,1) 0 meanEVREV.run.pyr.Rnomod(:,2)],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');
    xlabel('RUN Pyr Only Rip Mod/NoMod');
    ylabel(['n=' int2str(length(sess.run.pyr.Rmod)) 'sess']);
    sigstar({[1 3],[2 4]},[p.run.pyr.Rmod,p.run.pyr.Rnomod]);
  end
  subplot(2,2,3);
  barwitherr(semEVREV.run.pyr.Rall(:,1:2),meanEVREV.run.pyr.Rall(:,1:2),'FaceColor',rgb('Crimson'),'EdgeColor','none');
  set(gca,'XTickLabel',titrun);
  xlabel('RUN Pyr Only');
  ylabel(['n=' int2str(length(sess.run.pyr.Rall)) 'sess']);
  sigstar([1 2],p.run.pyr.Rall);
  suptitle(suptit);
end

if strcmp(EVtype,'cross')
  evtypetitle='Cross-Structure';
elseif strcmp(EVtype,'intra')
  evtypetitle='Intra-Structure';
end

savefig=input('Save Figures?','s')
if strcmp(savefig,'yes')
  cd('/media/Data-01/All-Rats/AllRats-ExplainedVariance')
  figure(f1)
  plot2svg(['Figure1-' struct '-' evtypetitle '-SWS.svg'],gcf);
  saveas(gcf,['Figure1-' struct '-' evtypetitle '-SWS'],'png');
  saveas(gcf,['Figure1-' struct '-' evtypetitle '-SWS'],'pdf');
  figure(f2)
  plot2svg(['Figure2-' struct '-' evtypetitle '-REM.svg'],gcf);
  saveas(gcf,['Figure2-' struct '-' evtypetitle '-REM'],'png');
  saveas(gcf,['Figure2-' struct '-' evtypetitle '-REM'],'pdf');
  figure(f3)
  plot2svg(['Figure3-' struct '-' evtypetitle '-RIP.svg'],gcf);
  saveas(gcf,['Figure3-' struct '-' evtypetitle '-RIP'],'png');
  saveas(gcf,['Figure3-' struct '-' evtypetitle '-RIP'],'pdf');
%    figure(f4)
%    plot2svg(['Figure4-' struct '-' evtypetitle '-RUN.svg'],gcf);
%    saveas(gcf,['Figure4-' struct '-' evtypetitle '-RUN'],'png');
%    saveas(gcf,['Figure4-' struct '-' evtypetitle '-RUN'],'pdf');
end
if strcmp(savevar,'on')
  cd('/media/Data-01/All-Rats/AllRats-ExplainedVariance')
  save(['EVREV-Ultimate-' struct '.mat'],'EVREV','sess','p','stats');
end

%  figure;
%  EVratio.run=(EVREV.run.pyr.Rall(:,1)-EVREV.run.pyr.Rall(:,2))./(EVREV.run.pyr.Rall(:,1)+EVREV.run.pyr.Rall(:,2));
%  EVratio.sws=(EVREV.sws.pyr.Rall(2:end,1)-EVREV.sws.pyr.Rall(2:end,2))./(EVREV.sws.pyr.Rall(2:end,1)+EVREV.sws.pyr.Rall(2:end,2));
%  plot(EVratio.run,EVratio.sws,'.')
%  lsline;
%  
%  EVgain.run=EVREV.run.pyr.Rall(:,1)./EVREV.run.pyr.Rall(:,2);
%  EVgain.sws=EVREV.sws.pyr.Rall(2:end,1)./EVREV.sws.pyr.Rall(2:end,2);
%  figure;
%  plot(EVgain.run,EVgain.sws,'.')
%  lsline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  figure;
%  plot(allcorr.pre,allcorr.post,'.');
%  xlabel('PRE');ylabel('POST');
%  
%  ylim([-0.15 0.25]);
%  xlim([-0.15 0.25]);
%  hold on;
%  plot([-0.15 0.25],[-0.15 0.25],'k:');
%  [r1,p1]=corrcoef(allcorr.pre,allcorr.post,'rows','complete');
%  r1=r1(1,2);p1=p1(1,2);
%  pf=polyfit(allcorr.pre(~isnan(allcorr.pre)&~isnan(allcorr.post)),allcorr.post(~isnan(allcorr.pre)&~isnan(allcorr.post)),1);
%  plot([-0.15 0.25],[-0.15*pf(1)+pf(2) 0.25*pf(1)+pf(2)],'r');
%  suptitle([brainstate ' r:' num2str(r1) ' p:' num2str(p1)]);
%  
%  
%  figure;
%  subplot(1,3,2);hold on;
%  [h,bins]=hist(allcorr.run,100);
%  h2=hist(allcorr.run(allp.run<0.01),bins);
%  bar(bins,h,'LineStyle','none');
%  bar(bins,h2,'r','LineStyle','none')
%  xlabel(int2str(length(allcorr.run)));
%  xlim([-0.1 0.4])
%  subplot(1,3,1);hold on;
%  h=hist(allcorr.pre,bins);
%  h2=hist(allcorr.pre(allp.pre<0.01),bins);
%  bar(bins,h,'LineStyle','none');
%  bar(bins,h2,'r','LineStyle','none');
%  xlabel(int2str(length(allcorr.pre)));
%  xlim([-0.1 0.4])
%  subplot(1,3,3);hold on;
%  h=hist(allcorr.post,bins);
%  h2=hist(allcorr.post(allp.post<0.01),bins);
%  bar(bins,h,'LineStyle','none');
%  bar(bins,h2,'r','LineStyle','none');
%  xlabel(int2str(length(allcorr.post)));
%  xlim([-0.1 0.4])
%  suptitle(brainstate);
%  if strcmp(brainstate,'ripples')
%    figure;
%    subplot(1,3,2);hold on;
%    [h,bins]=hist(allcorr.run,100);
%    h2=hist(allcorr.run(allp.run<0.01),bins);
%    bar(bins,h,'LineStyle','none');
%    bar(bins,h2,'r','LineStyle','none');
%    xlabel(int2str(length(allcorr.run)));
%    xlim([-0.1 0.4])
%    subplot(1,3,1);hold on;
%    h=hist(allcorr.preout,bins);
%    h2=hist(allcorr.preout(allp.preout<0.01),bins);
%    bar(bins,h,'LineStyle','none');
%    bar(bins,h2,'r','LineStyle','none');
%    xlabel(int2str(length(allcorr.preout)));
%    xlim([-0.1 0.4])
%    subplot(1,3,3);hold on;
%    h=hist(allcorr.postout,bins);
%    h2=hist(allcorr.postout(allp.postout<0.01),bins);
%    bar(bins,h,'LineStyle','none');
%    bar(bins,h2,'r','LineStyle','none');
%    xlabel(int2str(length(allcorr.post)));
%    xlim([-0.1 0.4])
%    suptitle('Outside ripples');
%  end
%  
%  EVratio=(EVREV(:,1)-EVREV(:,2))./(EVREV(:,1)+EVREV(:,2));
%  EVgain=EVREV(:,1)./EVREV(:,2);
%  if strcmp(brainstate,'sws');
%    begEVratio=(EVREV(:,3)-EVREV(:,4))./(EVREV(:,3)+EVREV(:,4));
%    posEVratio=(EVREV(EVREV(:,3)-EVREV(:,4)>0,3)-EVREV(EVREV(:,3)-EVREV(:,4)>0,4))./(EVREV(EVREV(:,3)-EVREV(:,4)>0,3)+EVREV(EVREV(:,3)-EVREV(:,4)>0,4))
%  end
%  
%  figure('Position',[362 69 1402 894]);				% BEHAV : ratio.run ratio.postrun perf
%  subplot(3,4,1)
%  plot(EVratio(:,1),behav(:,1),'.');
%  [a,p]=corrcoef(EVratio,behav(:,1),'rows','complete');
%  xlabel('EVratio');
%  ylabel('Velocity ratio RUN');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  subplot(3,4,2)
%  plot(EVgain(:,1),behav(:,1),'.');
%  [a,p]=corrcoef(EVgain,behav(:,1),'rows','complete');
%  xlabel('EVgain');
%  ylabel('Velocity ratio RUN');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  subplot(3,4,5)
%  plot(EVratio(:,1),behav(:,2),'.');
%  [a,p]=corrcoef(EVratio,behav(:,2),'rows','complete');
%  xlabel('EVratio');
%  ylabel('Velocity ratio postRUN');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  subplot(3,4,6)
%  plot(EVgain(:,1),behav(:,2),'.');
%  [a,p]=corrcoef(EVgain,behav(:,2),'rows','complete');
%  xlabel('EVgain');
%  ylabel('Velocity ratio postRUN');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  subplot(3,4,9)
%  plot(EVratio(:,1),behav(:,3),'.');
%  [a,p]=corrcoef(EVratio,behav(:,3),'rows','complete');
%  xlabel('EVratio');
%  ylabel('Perf');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  subplot(3,4,10)
%  plot(EVgain(:,1),behav(:,3),'.');
%  [a,p]=corrcoef(EVgain,behav(:,3),'rows','complete');
%  xlabel('EVgain');
%  ylabel('Perf');
%  title([brainstate ' p=' num2str(p(1,2))]);
%  if strcmp(brainstate,'sws');
%    subplot(3,4,3)
%    plot(begEVratio(:,1),behav(:,1),'.');
%    [a,p]=corrcoef(begEVratio,behav(:,1),'rows','complete');
%    xlabel('begEVratio');
%    ylabel('Velocity ratio RUN');
%    title([brainstate ' p=' num2str(p(1,2))]);	
%    subplot(3,4,4)
%    plot(posEVratio(:,1),behav(EVREV(:,3)-EVREV(:,4)>0,1),'.');
%    [a,p]=corrcoef(posEVratio,behav(EVREV(:,3)-EVREV(:,4)>0,1),'rows','complete');
%    xlabel('posEVratio');
%    ylabel('Velocity ratio RUN');
%    title([brainstate ' p=' num2str(p(1,2))]);	
%    subplot(3,4,7)
%    plot(begEVratio(:,1),behav(:,2),'.');
%    [a,p]=corrcoef(begEVratio,behav(:,2),'rows','complete');
%    xlabel('begEVratio');
%    ylabel('Velocity ratio postRUN');
%    title([brainstate ' p=' num2str(p(1,2))]);	
%    subplot(3,4,8)
%    plot(posEVratio(:,1),behav(EVREV(:,3)-EVREV(:,4)>0,2),'.');
%    [a,p]=corrcoef(posEVratio,behav(EVREV(:,3)-EVREV(:,4)>0,2),'rows','complete');
%    xlabel('posEVratio');
%    ylabel('Velocity ratio postRUN');
%    title([brainstate ' p=' num2str(p(1,2))]);	
%    subplot(3,4,11)
%    plot(begEVratio(:,1),behav(:,3),'.');
%    [a,p]=corrcoef(begEVratio,behav(:,3),'rows','complete');
%    xlabel('begEVratio');
%    ylabel('Perf');
%    title([brainstate ' p=' num2str(p(1,2))]);	
%    subplot(3,4,12)
%    plot(posEVratio(:,1),behav(EVREV(:,3)-EVREV(:,4)>0,3),'.');
%    [a,p]=corrcoef(posEVratio,behav(EVREV(:,3)-EVREV(:,4)>0,3),'rows','complete');
%    xlabel('posEVratio');
%    ylabel('Perf');
%    title([brainstate ' p=' num2str(p(1,2))]);
%  end


