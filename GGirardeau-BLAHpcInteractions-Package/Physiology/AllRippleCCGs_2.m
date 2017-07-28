function allccg = AllRippleCCGs_2(session,presleep,postsleep)

%AllRippleCCGs_2 -Calculates spikes-ripples CCGs for all cells of a given session.
%
%  USAGE
%
%    allccg = AllRippleCCGs_2 (session, presleep, postsleep)
%
%    session            session path              
%    presleep           presleep epoch number
%    postsleep          postsleep epoch number
%
%  OUTPUT
%
%    allccg             matricx of cross-correolgrams between spikes and ripples (1 line/cell) for all cells of the session   
%
%  NOTE
%
%  SEE
%
%    See also : AllRippleCCGsAll
%
% Gabrielle Girardeau, 2016
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


cd (session);
xml=session(end-13:end);
SetCurrentSession(xml);

binSize1=0.001; % in seconds
duration1=0.4;
binSize2=0.01;
duration2=4;

load([xml '-States.mat']);
allpre=RunIntervals(presleep);
allpost=RunIntervals(postsleep);

Sws.all=sws;
Sws.pre=Restrict(sws,allpre);
Sws.post=Restrict(sws,allpost);

% initialize ccgvariables
allccg.ccg1=[];
allccg.gain1=[];
allccg.ccg2=[];
allccg.gain2=[];
allccg.t1=[];
allccg.t2=[];
allccg.binSize1=binSize1;
allccg.binSize2=binSize2;
allccg.duration1=duration1;
allccg.duration2=duration2;

allccg.pre.ccg1=[];
allccg.pre.gain1=[];
allccg.post.ccg2=[];
allccg.post.gain2=[];

% GetSpikes and compute index (spikes ; time/shank/cell/ID)
spikes1=GetSpikeTimes('output','numbered');
spikes2=GetSpikeTimes('output','full');
spikes=[spikes2 spikes1(:,2)];
spikes(spikes(:,3)==0,:)=[];
spikes(spikes(:,3)==1,:)=[];
% Remnumber IDs after removing MUA and artifacts
ids=unique(spikes(:,4));
for ii=1:length(ids)
  spikes(spikes(:,4)==ids(ii),4)=ii;
end
idx=unique(spikes(:,2:4),'rows'); % shank / cell/ ID

% Get Ripples and restrict to SWS
rip=GetRippleEvents;
ripples.all=Restrict(rip,Sws.all);
ripples.pre=Restrict(rip,Sws.pre);
ripples.post=Restrict(rip,Sws.post);

% Baseline epochs for baseline FR
bufferedripples=[rip(:,1)-0.1 rip(:,3)+0.1];
[baseline.all,ind]=SubtractIntervals(sws,bufferedripples,'strict','on'); % baseline = baseline intervals
baseline.pre=Restrict(baseline.all,allpre);
baseline.post=Restrict(baseline.all,allpost);

totalbaselinetime.all=sum(baseline.all(:,2)-baseline.all(:,1));
totalbaselinetime.pre=sum(baseline.pre(:,2)-baseline.pre(:,1));
totalbaselinetime.post=sum(baseline.post(:,2)-baseline.post(:,1));

baselinespikes.all=Restrict(spikes,baseline.all);
baselinespikes.pre=Restrict(spikes,baseline.pre);
baselinespikes.post=Restrict(spikes,baseline.post);

%baselineFR vector for all cells
for i=1:length(idx)
  baselineFR.all(i)=sum(baselinespikes.all(:,4)==i)./totalbaselinetime.all;
  baselineFR.pre(i)=sum(baselinespikes.pre(:,4)==i)./totalbaselinetime.pre;
  baselineFR.post(i)=sum(baselinespikes.post(:,4)==i)./totalbaselinetime.post;
end
baselineFR.all=baselineFR.all';
baselineFR.pre=baselineFR.pre';
baselineFR.post=baselineFR.post';

%prep for CCG
rippleccg.all=[ripples.all(:,2) ones(size(ripples.all,1),1)];
rippleccg.pre=[ripples.pre(:,2) ones(size(ripples.pre,1),1)];
rippleccg.post=[ripples.post(:,2) ones(size(ripples.post,1),1)];

spikesccg=spikes(:,[1 4]);

[s,ids,groups]=CCGParameters(rippleccg.all,1,spikesccg,2);
[ccg1,allccg.t1]=CCG(s,ids,'groups',groups,'binSize',binSize1,'duration',duration1,'smooth',1);
ccg1=squeeze(ccg1);
allccg.ccg1=ccg1';
allccg.ccg1=allccg.ccg1./(binSize1*size(ripples.all,1));
%
[ccg2,allccg.t2]=CCG(s,ids,'groups',groups,'binSize',binSize2,'duration',duration2,'smooth',1);
ccg2=squeeze(ccg2);
allccg.ccg2=ccg2';
allccg.ccg2=allccg.ccg2./(binSize2*size(ripples.all,1));
%%%
[s,ids,groups]=CCGParameters(rippleccg.pre,1,spikesccg,2);
[ccg1,allccg.t1]=CCG(s,ids,'groups',groups,'binSize',binSize1,'duration',duration1,'smooth',1);
ccg1=squeeze(ccg1);
allccg.pre.ccg1=ccg1';
allccg.pre.ccg1=allccg.pre.ccg1./(binSize1*size(ripples.pre,1));
%
[ccg2,allccg.t2]=CCG(s,ids,'groups',groups,'binSize',binSize2,'duration',duration2,'smooth',1);
ccg2=squeeze(ccg2);
allccg.pre.ccg2=ccg2';
allccg.pre.ccg2=allccg.pre.ccg2./(binSize2*size(ripples.pre,1));
%%%
[s,ids,groups]=CCGParameters(rippleccg.post,1,spikesccg,2);
[ccg1,allccg.t1]=CCG(s,ids,'groups',groups,'binSize',binSize1,'duration',duration1,'smooth',1);
ccg1=squeeze(ccg1);
allccg.post.ccg1=ccg1';
allccg.post.ccg1=allccg.post.ccg1./(binSize1*size(ripples.post,1));
%
[ccg2,allccg.t2]=CCG(s,ids,'groups',groups,'binSize',binSize2,'duration',duration2,'smooth',1);
ccg2=squeeze(ccg2);
allccg.post.ccg2=ccg2';
allccg.post.ccg2=allccg.post.ccg2./(binSize2*size(ripples.post,1));

matrixFR1.all=repmat(baselineFR.all,[1 size(allccg.ccg1,2)]);
matrixFR2.all=repmat(baselineFR.all,[1 size(allccg.ccg2,2)]);
matrixFR1.pre=repmat(baselineFR.pre,[1 size(allccg.pre.ccg1,2)]);
matrixFR2.pre=repmat(baselineFR.pre,[1 size(allccg.pre.ccg2,2)]);
matrixFR1.post=repmat(baselineFR.post,[1 size(allccg.post.ccg1,2)]);
matrixFR2.post=repmat(baselineFR.post,[1 size(allccg.post.ccg2,2)]);

allccg.gain1=allccg.ccg1./matrixFR1.all;
allccg.gain2=allccg.ccg2./matrixFR2.all;
allccg.pre.gain1=allccg.pre.ccg1./matrixFR1.pre;
allccg.pre.gain2=allccg.pre.ccg2./matrixFR2.pre;
allccg.post.gain1=allccg.post.ccg1./matrixFR1.post;
allccg.post.gain2=allccg.post.ccg2./matrixFR2.post;

imagesc(allccg.gain1);
figure;
imagesc(allccg.gain2);

save([xml '-RippleCCGs'],'idx','allccg');

%    % loop on cells
%    for j=1:size(cells,1)
%      cellref=[sigsessions(i,:) cells(j,:)];	% complete cellref rat/session/shank/number to retriece pInc/pDec
%      spikes=GetSpikeTimes(cells(j,:));
%      baselinespikes=Restrict(spikes,baseline);
%      baselinerate=length(baselinespikes)/sum(baseline(:,2)-baseline(:,1));
%      % spikecounts per bin, peri-ripple
%      spikesccg=[spikes 2*ones(length(spikes),1)];
%      forccg=[rippleccg; spikesccg];
%      forccg=sortrows(forccg,1);
%      % CCG short timescale
%      [ccg1,ti1]=CCG(forccg(:,1),forccg(:,2),'binSize',binSize1,'duration',0.4,'smooth',1);
%      ccg1=ccg1(:,1,2);
%      instccg1=ccg1./(binSize1*size(ripples,1));
%      gainccg1=instccg1./baselinerate;
%
%      normgainccg1=ZeroToOne(gainccg1);
%
%      % CCG long timescale
%  %      [sum2,tih2]=SyncHist(syncspk,indices,'durations',[-2 2],'nBins',401,'smooth',1);
%  %      ccg2=sum2;
%      [ccg2,ti2]=CCG(forccg(:,1),forccg(:,2),'binSize',binSize2,'duration',4,'smooth',1);
%      ccg2=ccg2(:,1,2);
%
%      instccg2=ccg2./(binSize2*size(ripples,1));
%      gainccg2=instccg2./baselinerate;
%
%      normgainccg2=ZeroToOne(gainccg2);
%
%  %      figure; bar(ti2,Smooth(gainccg2,1));
%      if sigcells(ismember(sigcells(:,1:4),cellref,'rows'),6)<pval
%        allccg.Inc1=[allccg.Inc1;gainccg1'];
%        allccg.Inc2=[allccg.Inc2;gainccg2'];
%        allccg.normInc1=[allccg.normInc1;normgainccg1'];
%        allccg.normInc2=[allccg.normInc1;normgainccg2'];
%        allccg.Inccells=[allccg.Inccells;sigcells(ismember(sigcells(:,1:4),cellref,'rows'),:)];
%      else
%        allccg.Dec1=[allccg.Dec1;gainccg1'];
%        allccg.Dec2=[allccg.Dec2;gainccg2'];
%        allccg.normDec1=[allccg.normDec1;normgainccg1'];
%        allccg.normDec2=[allccg.normDec1;normgainccg2'];
%        allccg.Deccells=[allccg.Deccells;sigcells(ismember(sigcells(:,1:4),cellref,'rows'),:)];
%      end
%    end
%    allccg.t1=ti1;
%    allccg.t2=ti2;
%  end
%
%  figure;
%  imagesc(allccg.Inc1);
%  xlabel(['Increasing cells - binSize: ' num2str(binSize1) ' duration: +/-' num2str(duration1/2) ' sec']);
%  figure;
%  imagesc(allccg.Inc2);
%  xlabel(['Increasing cells - binSize: ' num2str(binSize2) ' duration: +/-' num2str(duration2/2) ' sec']);
%  figure;
%  imagesc(allccg.Dec1);
%  xlabel(['Decreasing cells - binSize: ' num2str(binSize1) ' duration: +/-' num2str(duration1/2) ' sec']);
%  figure;
%  imagesc(allccg.Dec2);
%  xlabel(['Decreasing cells - binSize: ' num2str(binSize2) ' duration: +/-' num2str(duration2/2) ' sec']);
%
%  figure;
%  imagesc(allccg.normInc1);
%  xlabel(['Increasing cells - binSize: ' num2str(binSize1) ' duration: +/-' num2str(duration1/2) ' sec']);
%  figure;
%  imagesc(allccg.normInc2);
%  xlabel(['Increasing cells - binSize: ' num2str(binSize2) ' duration: +/-' num2str(duration2/2) ' sec']);
%  figure;
%  imagesc(allccg.normDec1);
%  xlabel(['Decreasing cells - binSize: ' num2str(binSize1) ' duration: +/-' num2str(duration1/2) ' sec']);
%  figure;
%  imagesc(allccg.normDec2);
%  xlabel(['Decreasing cells - binSize: ' num2str(binSize2) ' duration: +/-' num2str(duration2/2) ' sec']);
%
