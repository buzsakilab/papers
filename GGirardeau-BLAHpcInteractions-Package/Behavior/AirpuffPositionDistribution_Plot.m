function AirpuffPositionDistribution_Plot

%  AirpuffPositionDistribution_Plot - Plots behavioral measures against normalized airpuff position on the track 
%
%  USAGE
%
%    AirpuffPositionDistribution_Plot
%
%  NOTE
%
%    Adapted to specific data structure and session subselection. Uses stored variable "BehavioralMeasures.mat" and "NormAPpos.mat"
%
%  OUTPUT
%
%    Various plots
%
%  NOTE
%  
%  SEE
%
%    See also : AirpuffPositionDistribution, ZeroToOne (FMAToolbox)
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


cd('/media/Data-01/All-Rats/BehavioralMeasures');
load('BehavioralMeasures.mat');
load('NormAPpos.mat');

fullsess=ratsession(ismember(ratsession,ratsess,'rows'),:);
is.rat55=ismember(ratsession,fullsess,'rows');
is.rat53=ismember(ratsess,fullsess,'rows');

for i=1:length(PrePmeanVAll)
    RunVel(i)=PrePmeanVAll{i}.run;
    PrerunVel(i)=PrePmeanVAll{i}.prerun;
    PostrunVel(i)=PrePmeanVAll{i}.postrun;
end


RunVel=RunVel(is.rat55)';
PrerunVel=PrerunVel(is.rat55)';
PostrunVel=PostrunVel(is.rat55)';

NormAP=NormalizedAPpos(is.rat53);
dirs=directions(is.rat53);

flippedNormAP=NormAP;
flippedNormAP(dirs==1)=abs(flippedNormAP(dirs==1)-1)

%  figure;
%  subplot(2,2,1);hold on;
%  plot(NormAP(fullsess(:,1)==8),RunVel(fullsess(:,1)==8),'r.');
%  plot(NormAP(fullsess(:,1)==8),PrerunVel(fullsess(:,1)==8),'g.');
%  plot(NormAP(fullsess(:,1)==8),PostrunVel(fullsess(:,1)==8),'b.');
%  ylim([0 90])
%  xlim([0 1])
%  subplot(2,2,2);hold on;
%  plot(NormAP(fullsess(:,1)==9),RunVel(fullsess(:,1)==9),'r.');
%  plot(NormAP(fullsess(:,1)==9),PrerunVel(fullsess(:,1)==9),'g.');
%  plot(NormAP(fullsess(:,1)==9),PostrunVel(fullsess(:,1)==9),'b.');
%  ylim([0 90])
%  xlim([0 1])
%  subplot(2,2,3);hold on;
%  plot(NormAP(fullsess(:,1)==10),RunVel(fullsess(:,1)==10),'r.');
%  plot(NormAP(fullsess(:,1)==10),PrerunVel(fullsess(:,1)==10),'g.');
%  plot(NormAP(fullsess(:,1)==10),PostrunVel(fullsess(:,1)==10),'b.');
%  ylim([0 90])
%  xlim([0 1])
%  subplot(2,2,4);hold on;
%  xlim([0 1])
%  plot(NormAP(fullsess(:,1)==11),RunVel(fullsess(:,1)==11),'r.');
%  plot(NormAP(fullsess(:,1)==11),PrerunVel(fullsess(:,1)==11),'g.');
%  plot(NormAP(fullsess(:,1)==11),PostrunVel(fullsess(:,1)==11),'b.');
%  ylim([0 90])
%  xlim([0 1])

figure;hold on;
plot(NormAP,RunVel,'r.');
plot(NormAP,PrerunVel,'g.');
plot(NormAP,PostrunVel,'b.');

figure;
subplot(1,3,1);
plot(flippedNormAP,PrerunVel,'.b');
lsline
ylim([0 100])
subplot(1,3,2);
plot(flippedNormAP,RunVel,'.k');
lsline
ylim([0 100])
subplot(1,3,3);
plot(flippedNormAP,PostrunVel,'.r');
lsline
ylim([0 100])

figure;
plot([flippedNormAP;flippedNormAP],[PrerunVel;PostrunVel],'.b');
lsline

[r1,p1]=corrcoef(flippedNormAP,PrerunVel)
[r2,p2]=corrcoef(flippedNormAP,RunVel)
[r3,p3]=corrcoef(flippedNormAP(~isnan(PostrunVel)),PostrunVel(~isnan(PostrunVel)))


figure;
subplot(1,3,1);
plot(NormAP,PrerunVel,'.b');
lsline
ylim([0 100])
subplot(1,3,2);
plot(NormAP,RunVel,'.k');
lsline
ylim([0 100])
subplot(1,3,3);
plot(NormAP,PostrunVel,'.r');
lsline
ylim([0 100])

CenterDist=abs(NormAP-0.5)

figure;
subplot(1,3,1);
plot(CenterDist,PrerunVel,'.b');
lsline
ylim([0 100])
subplot(1,3,2);
plot(CenterDist,RunVel,'.k');
lsline
ylim([0 100])
subplot(1,3,3);
plot(CenterDist,PostrunVel,'.r');
lsline
ylim([0 100])

[r1,p1]=corrcoef(CenterDist,PrerunVel)
[r2,p2]=corrcoef(CenterDist,RunVel)
[r3,p3]=corrcoef(CenterDist(~isnan(PostrunVel)),PostrunVel(~isnan(PostrunVel)))

QuartNormAP=NormAP;
QuartNormAP=QuartNormAP*10;
QuartNormAP=floor(QuartNormAP);

figure;hold on;
subplot(1,3,2);hold on;
plot(QuartNormAP,RunVel,'r.');
boxplot(RunVel,QuartNormAP);
xlim([0 9])
ylim([0 110])
subplot(1,3,1);hold on;
plot(QuartNormAP,PrerunVel,'g.');
boxplot(PrerunVel,QuartNormAP);
xlim([0 9])
ylim([0 110])
subplot(1,3,3);hold on;
plot(QuartNormAP,PostrunVel,'b.');
boxplot(PostrunVel,QuartNormAP);
xlim([0 9])
ylim([0 110])


FourNormAP=QuartNormAP;
FourNormAP(FourNormAP==5)=4;
FourNormAP(FourNormAP==6)=3;
FourNormAP(FourNormAP==7)=2;
FourNormAP(FourNormAP==8)=1;

figure;hold on;
subplot(1,3,2);hold on;
plot(FourNormAP,RunVel,'r.');
boxplot(RunVel,FourNormAP);
xlim([0 5])
ylim([0 110])
subplot(1,3,1);hold on;
plot(FourNormAP,PrerunVel,'g.');
boxplot(PrerunVel,FourNormAP);
xlim([0 5])
ylim([0 110])
subplot(1,3,3);hold on;
plot(FourNormAP,PostrunVel,'b.');
boxplot(PostrunVel,FourNormAP);
xlim([0 5])
ylim([0 110])

[p,tab,stats]=kruskalwallis(RunVel,FourNormAP);
[c,m,h]=multcompare(stats);
