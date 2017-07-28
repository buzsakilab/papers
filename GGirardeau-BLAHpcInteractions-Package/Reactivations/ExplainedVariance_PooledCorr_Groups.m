function [EV,REV] = ExplainedVariance_PooledCorr_Groups

% Explained Variance_PooledCorr_Groups - Calculates and plots EV and REVs for different groups of Hpc-BLA pairs (based on ripple-modulation of BLA cell, RUN correlation and firing rates) 
%  
%  USAGE
%
%    [EV, REV] = ExplainedVariance_PooledCorr_Groups
%
%    struct		structure name (ex : 'BLA')
%    threshold		threshold for contributing and anticontributing cells
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
%  NOTE
%
%  SEE
%
%    See also : binspikes 
%
% March 2016, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/AllRats-EVcontribs.mat')

bothpyr=sum(Index.Type,2)==2;

ripple.up=Index.RippleMod==1;
ripple.down=Index.RippleMod==2;
ripple.none=Index.RippleMod==0;

run.sigup=Allpval.run(:,2)<0.01&Allcorr.run(:,2)>0;
run.sigdown=Allpval.run(:,2)<0.01&Allcorr.run(:,2)<0;
run.signot=Allpval.run(:,2)>=0.01&Allcorr.run(:,2)<0;
% tests for run significance automatically removes NaN.

%%%%
ripup.runup=bothpyr&ripple.up&run.sigup;
ripup.rundown=bothpyr&ripple.up&run.sigdown;
ripup.runnot=bothpyr&ripple.up&run.signot;
ripup.all=bothpyr&ripple.up;

ripdown.runup=bothpyr&ripple.down&run.sigup;
ripdown.rundown=bothpyr&ripple.down&run.sigdown;
ripdown.runnot=bothpyr&ripple.down&run.signot;
ripdown.all=bothpyr&ripple.down;

ripnone.runup=bothpyr&ripple.none&run.sigup;
ripnone.rundown=bothpyr&ripple.none&run.sigdown;
ripnone.runnot=bothpyr&ripple.none&run.signot;
ripnone.all=bothpyr&ripple.none;

% Need to remove NaNs for ripple-mod only EV/REV 
nonan=~sum(isnan(Allcorr.sws),2)>0;

%%%% Ripple up
% Run sig pos
r_runpost.ripup.runup=corrcoef(Allcorr.sws(ripup.runup,3),Allcorr.sws(ripup.runup,2));r_runpost.ripup.runup=r_runpost.ripup.runup(1,2);
r_runpre.ripup.runup=corrcoef(Allcorr.sws(ripup.runup,1),Allcorr.sws(ripup.runup,2));r_runpre.ripup.runup=r_runpre.ripup.runup(1,2);
r_prepost.ripup.runup=corrcoef(Allcorr.sws(ripup.runup,1),Allcorr.sws(ripup.runup,3));r_prepost.ripup.runup=r_prepost.ripup.runup(1,2);
EV.ripup.runup=((r_runpost.ripup.runup-r_runpre.ripup.runup*r_prepost.ripup.runup)/sqrt((1-r_runpre.ripup.runup^2)*(1-r_prepost.ripup.runup^2)))^2;
REV.ripup.runup=((r_runpre.ripup.runup-r_runpost.ripup.runup*r_prepost.ripup.runup)/sqrt((1-r_runpost.ripup.runup^2)*(1-r_prepost.ripup.runup^2)))^2;

% Run sig neg
r_runpost.ripup.rundown=corrcoef(Allcorr.sws(ripup.rundown,3),Allcorr.sws(ripup.rundown,2));r_runpost.ripup.rundown=r_runpost.ripup.rundown(1,2);
r_runpre.ripup.rundown=corrcoef(Allcorr.sws(ripup.rundown,1),Allcorr.sws(ripup.rundown,2));r_runpre.ripup.rundown=r_runpre.ripup.rundown(1,2);
r_prepost.ripup.rundown=corrcoef(Allcorr.sws(ripup.rundown,1),Allcorr.sws(ripup.rundown,3));r_prepost.ripup.rundown=r_prepost.ripup.rundown(1,2);
EV.ripup.rundown=((r_runpost.ripup.rundown-r_runpre.ripup.rundown*r_prepost.ripup.rundown)/sqrt((1-r_runpre.ripup.rundown^2)*(1-r_prepost.ripup.rundown^2)))^2;
REV.ripup.rundown=((r_runpre.ripup.rundown-r_runpost.ripup.rundown*r_prepost.ripup.rundown)/sqrt((1-r_runpost.ripup.rundown^2)*(1-r_prepost.ripup.rundown^2)))^2;

% Run none
r_runpost.ripup.runnot=corrcoef(Allcorr.sws(ripup.runnot,3),Allcorr.sws(ripup.runnot,2));r_runpost.ripup.runnot=r_runpost.ripup.runnot(1,2);
r_runpre.ripup.runnot=corrcoef(Allcorr.sws(ripup.runnot,1),Allcorr.sws(ripup.runnot,2));r_runpre.ripup.runnot=r_runpre.ripup.runnot(1,2);
r_prepost.ripup.runnot=corrcoef(Allcorr.sws(ripup.runnot,1),Allcorr.sws(ripup.runnot,3));r_prepost.ripup.runnot=r_prepost.ripup.runnot(1,2);
EV.ripup.runnot=((r_runpost.ripup.runnot-r_runpre.ripup.runnot*r_prepost.ripup.runnot)/sqrt((1-r_runpre.ripup.runnot^2)*(1-r_prepost.ripup.runnot^2)))^2;
REV.ripup.runnot=((r_runpre.ripup.runnot-r_runpost.ripup.runnot*r_prepost.ripup.runnot)/sqrt((1-r_runpost.ripup.runnot^2)*(1-r_prepost.ripup.runnot^2)))^2;

% All ripple up
r_runpost.ripup.all=corrcoef(Allcorr.sws(ripup.all&nonan,3),Allcorr.sws(ripup.all&nonan,2));r_runpost.ripup.all=r_runpost.ripup.all(1,2);
r_runpre.ripup.all=corrcoef(Allcorr.sws(ripup.all&nonan,1),Allcorr.sws(ripup.all&nonan,2));r_runpre.ripup.all=r_runpre.ripup.all(1,2);
r_prepost.ripup.all=corrcoef(Allcorr.sws(ripup.all&nonan,1),Allcorr.sws(ripup.all&nonan,3));r_prepost.ripup.all=r_prepost.ripup.all(1,2);
EV.ripup.all=((r_runpost.ripup.all-r_runpre.ripup.all*r_prepost.ripup.all)/sqrt((1-r_runpre.ripup.all^2)*(1-r_prepost.ripup.all^2)))^2;
REV.ripup.all=((r_runpre.ripup.all-r_runpost.ripup.all*r_prepost.ripup.all)/sqrt((1-r_runpost.ripup.all^2)*(1-r_prepost.ripup.all^2)))^2;


%%%% Ripple down
% Run sig pos
r_runpost.ripdown.runup=corrcoef(Allcorr.sws(ripdown.runup,3),Allcorr.sws(ripdown.runup,2));r_runpost.ripdown.runup=r_runpost.ripdown.runup(1,2);
r_runpre.ripdown.runup=corrcoef(Allcorr.sws(ripdown.runup,1),Allcorr.sws(ripdown.runup,2));r_runpre.ripdown.runup=r_runpre.ripdown.runup(1,2);
r_prepost.ripdown.runup=corrcoef(Allcorr.sws(ripdown.runup,1),Allcorr.sws(ripdown.runup,3));r_prepost.ripdown.runup=r_prepost.ripdown.runup(1,2);
EV.ripdown.runup=((r_runpost.ripdown.runup-r_runpre.ripdown.runup*r_prepost.ripdown.runup)/sqrt((1-r_runpre.ripdown.runup^2)*(1-r_prepost.ripdown.runup^2)))^2;
REV.ripdown.runup=((r_runpre.ripdown.runup-r_runpost.ripdown.runup*r_prepost.ripdown.runup)/sqrt((1-r_runpost.ripdown.runup^2)*(1-r_prepost.ripdown.runup^2)))^2;

% Run sig neg
r_runpost.ripdown.rundown=corrcoef(Allcorr.sws(ripdown.rundown,3),Allcorr.sws(ripdown.rundown,2));r_runpost.ripdown.rundown=r_runpost.ripdown.rundown(1,2);
r_runpre.ripdown.rundown=corrcoef(Allcorr.sws(ripdown.rundown,1),Allcorr.sws(ripdown.rundown,2));r_runpre.ripdown.rundown=r_runpre.ripdown.rundown(1,2);
r_prepost.ripdown.rundown=corrcoef(Allcorr.sws(ripdown.rundown,1),Allcorr.sws(ripdown.rundown,3));r_prepost.ripdown.rundown=r_prepost.ripdown.rundown(1,2);
EV.ripdown.rundown=((r_runpost.ripdown.rundown-r_runpre.ripdown.rundown*r_prepost.ripdown.rundown)/sqrt((1-r_runpre.ripdown.rundown^2)*(1-r_prepost.ripdown.rundown^2)))^2;
REV.ripdown.rundown=((r_runpre.ripdown.rundown-r_runpost.ripdown.rundown*r_prepost.ripdown.rundown)/sqrt((1-r_runpost.ripdown.rundown^2)*(1-r_prepost.ripdown.rundown^2)))^2;

% Run none
r_runpost.ripdown.runnot=corrcoef(Allcorr.sws(ripdown.runnot,3),Allcorr.sws(ripdown.runnot,2));r_runpost.ripdown.runnot=r_runpost.ripdown.runnot(1,2);
r_runpre.ripdown.runnot=corrcoef(Allcorr.sws(ripdown.runnot,1),Allcorr.sws(ripdown.runnot,2));r_runpre.ripdown.runnot=r_runpre.ripdown.runnot(1,2);
r_prepost.ripdown.runnot=corrcoef(Allcorr.sws(ripdown.runnot,1),Allcorr.sws(ripdown.runnot,3));r_prepost.ripdown.runnot=r_prepost.ripdown.runnot(1,2);
EV.ripdown.runnot=((r_runpost.ripdown.runnot-r_runpre.ripdown.runnot*r_prepost.ripdown.runnot)/sqrt((1-r_runpre.ripdown.runnot^2)*(1-r_prepost.ripdown.runnot^2)))^2;
REV.ripdown.runnot=((r_runpre.ripdown.runnot-r_runpost.ripdown.runnot*r_prepost.ripdown.runnot)/sqrt((1-r_runpost.ripdown.runnot^2)*(1-r_prepost.ripdown.runnot^2)))^2;

% All Ripple down
r_runpost.ripdown.all=corrcoef(Allcorr.sws(ripdown.all&nonan,3),Allcorr.sws(ripdown.all&nonan,2));r_runpost.ripdown.all=r_runpost.ripdown.all(1,2);
r_runpre.ripdown.all=corrcoef(Allcorr.sws(ripdown.all&nonan,1),Allcorr.sws(ripdown.all&nonan,2));r_runpre.ripdown.all=r_runpre.ripdown.all(1,2);
r_prepost.ripdown.all=corrcoef(Allcorr.sws(ripdown.all&nonan,1),Allcorr.sws(ripdown.all&nonan,3));r_prepost.ripdown.all=r_prepost.ripdown.all(1,2);
EV.ripdown.all=((r_runpost.ripdown.all-r_runpre.ripdown.all*r_prepost.ripdown.all)/sqrt((1-r_runpre.ripdown.all^2)*(1-r_prepost.ripdown.all^2)))^2;
REV.ripdown.all=((r_runpre.ripdown.all-r_runpost.ripdown.all*r_prepost.ripdown.all)/sqrt((1-r_runpost.ripdown.all^2)*(1-r_prepost.ripdown.all^2)))^2;

%%%% Ripple none
% Run sig pos
r_runpost.ripnone.runup=corrcoef(Allcorr.sws(ripnone.runup,3),Allcorr.sws(ripnone.runup,2));r_runpost.ripnone.runup=r_runpost.ripnone.runup(1,2);
r_runpre.ripnone.runup=corrcoef(Allcorr.sws(ripnone.runup,1),Allcorr.sws(ripnone.runup,2));r_runpre.ripnone.runup=r_runpre.ripnone.runup(1,2);
r_prepost.ripnone.runup=corrcoef(Allcorr.sws(ripnone.runup,1),Allcorr.sws(ripnone.runup,3));r_prepost.ripnone.runup=r_prepost.ripnone.runup(1,2);
EV.ripnone.runup=((r_runpost.ripnone.runup-r_runpre.ripnone.runup*r_prepost.ripnone.runup)/sqrt((1-r_runpre.ripnone.runup^2)*(1-r_prepost.ripnone.runup^2)))^2;
REV.ripnone.runup=((r_runpre.ripnone.runup-r_runpost.ripnone.runup*r_prepost.ripnone.runup)/sqrt((1-r_runpost.ripnone.runup^2)*(1-r_prepost.ripnone.runup^2)))^2;

% Run sig neg
r_runpost.ripnone.rundown=corrcoef(Allcorr.sws(ripnone.rundown,3),Allcorr.sws(ripnone.rundown,2));r_runpost.ripnone.rundown=r_runpost.ripnone.rundown(1,2);
r_runpre.ripnone.rundown=corrcoef(Allcorr.sws(ripnone.rundown,1),Allcorr.sws(ripnone.rundown,2));r_runpre.ripnone.rundown=r_runpre.ripnone.rundown(1,2);
r_prepost.ripnone.rundown=corrcoef(Allcorr.sws(ripnone.rundown,1),Allcorr.sws(ripnone.rundown,3));r_prepost.ripnone.rundown=r_prepost.ripnone.rundown(1,2);
EV.ripnone.rundown=((r_runpost.ripnone.rundown-r_runpre.ripnone.rundown*r_prepost.ripnone.rundown)/sqrt((1-r_runpre.ripnone.rundown^2)*(1-r_prepost.ripnone.rundown^2)))^2;
REV.ripnone.rundown=((r_runpre.ripnone.rundown-r_runpost.ripnone.rundown*r_prepost.ripnone.rundown)/sqrt((1-r_runpost.ripnone.rundown^2)*(1-r_prepost.ripnone.rundown^2)))^2;

% Run none
r_runpost.ripnone.runnot=corrcoef(Allcorr.sws(ripnone.runnot,3),Allcorr.sws(ripnone.runnot,2));r_runpost.ripnone.runnot=r_runpost.ripnone.runnot(1,2);
r_runpre.ripnone.runnot=corrcoef(Allcorr.sws(ripnone.runnot,1),Allcorr.sws(ripnone.runnot,2));r_runpre.ripnone.runnot=r_runpre.ripnone.runnot(1,2);
r_prepost.ripnone.runnot=corrcoef(Allcorr.sws(ripnone.runnot,1),Allcorr.sws(ripnone.runnot,3));r_prepost.ripnone.runnot=r_prepost.ripnone.runnot(1,2);
EV.ripnone.runnot=((r_runpost.ripnone.runnot-r_runpre.ripnone.runnot*r_prepost.ripnone.runnot)/sqrt((1-r_runpre.ripnone.runnot^2)*(1-r_prepost.ripnone.runnot^2)))^2;
REV.ripnone.runnot=((r_runpre.ripnone.runnot-r_runpost.ripnone.runnot*r_prepost.ripnone.runnot)/sqrt((1-r_runpost.ripnone.runnot^2)*(1-r_prepost.ripnone.runnot^2)))^2;

% All ripple none
r_runpost.ripnone.all=corrcoef(Allcorr.sws(ripnone.all&nonan,3),Allcorr.sws(ripnone.all&nonan,2));r_runpost.ripnone.all=r_runpost.ripnone.all(1,2);
r_runpre.ripnone.all=corrcoef(Allcorr.sws(ripnone.all&nonan,1),Allcorr.sws(ripnone.all&nonan,2));r_runpre.ripnone.all=r_runpre.ripnone.all(1,2);
r_prepost.ripnone.all=corrcoef(Allcorr.sws(ripnone.all&nonan,1),Allcorr.sws(ripnone.all&nonan,3));r_prepost.ripnone.all=r_prepost.ripnone.all(1,2);
EV.ripnone.all=((r_runpost.ripnone.all-r_runpre.ripnone.all*r_prepost.ripnone.all)/sqrt((1-r_runpre.ripnone.all^2)*(1-r_prepost.ripnone.all^2)))^2;
REV.ripnone.all=((r_runpre.ripnone.all-r_runpost.ripnone.all*r_prepost.ripnone.all)/sqrt((1-r_runpost.ripnone.all^2)*(1-r_prepost.ripnone.all^2)))^2;


%% PLOT
figure;hold on;
bar([EV.ripup.all REV.ripup.all 0 0 0 0 ],'FaceColor',rgb('Crimson'),'EdgeColor','none');hold on;
bar([0 0 EV.ripnone.all REV.ripnone.all 0 0],'FaceColor',rgb('DarkGrey'),'EdgeColor','none');hold on;
bar([0 0 0 0 EV.ripdown.all REV.ripdown.all],'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');hold on;
xlabel('EV for ripple modulation groups, all pyr')


f1=figure;hold on;
bar([EV.ripup.runup REV.ripup.runup 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('Crimson'),'EdgeColor','none');hold on;
bar([0 0 EV.ripup.runnot REV.ripup.runnot 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkGrey'),'EdgeColor','none');hold on;
bar([0 0 0 0 EV.ripup.rundown REV.ripup.rundown 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('FireBrick'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 EV.ripnone.runup REV.ripnone.runup 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('SeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 EV.ripnone.runnot REV.ripnone.runnot 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 EV.ripnone.rundown REV.ripnone.rundown 0 0 0 0 0 0],'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 EV.ripdown.runup REV.ripdown.runup 0 0 0 0],'FaceColor',rgb('SteelBlue'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 0 0 EV.ripdown.runnot REV.ripdown.runnot 0 0],'FaceColor',rgb('LightBlue'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 EV.ripdown.rundown REV.ripdown.rundown],'FaceColor',rgb('DodgerBlue'),'EdgeColor','none');


MeanPairFR=mean(Index.FR,2);
PyrMeanFR=MeanPairFR(bothpyr);
MeanFRoct=quantile(PyrMeanFR,7);
MeanFRgroup=QuantileGroups(MeanPairFR,MeanFRoct);

% Octiles can be calculated either on FR pooled per pair or per cell. (per pair = Given firing rate appears the number of time that the cell has a pair partner.)
PyrBLAFR=Index.FR(bothpyr,2);
%  PyrBLAFR=unique(Index.FR(bothpyr,2));
BLAFRoct=quantile(PyrBLAFR,7);
BLAFRgroups=QuantileGroups(Index.FR(:,2),BLAFRoct);

PyrHpcFR=Index.FR(bothpyr,1);
%  PyrHpcFR=unique(Index.FR(bothpyr,1));
HpcFRoct=quantile(PyrHpcFR,7);
HpcFRgroups=QuantileGroups(Index.FR(:,1),HpcFRoct);


nonan=~isnan(Allcorr.sws(:,1))&~isnan(Allcorr.sws(:,2))&~isnan(Allcorr.sws(:,3));

for i=1:8
  R_runpost=corrcoef(Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,3),Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,2));r_runpost.rate(i)=R_runpost(1,2);
  R_runpre=corrcoef(Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,1),Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,2));r_runpre.rate(i)=R_runpre(1,2);
  R_prepost=corrcoef(Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,1),Allcorr.sws(bothpyr&BLAFRgroups==i&nonan,3));r_prepost.rate(i)=R_prepost(1,2);
  EV.rate(i)=((r_runpost.rate(i)-r_runpre.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpre.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
  REV.rate(i)=((r_runpre.rate(i)-r_runpost.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpost.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
end  
 
f2=figure('Position',[488 68 1004 935]);hold on;
subplot(3,2,1)
bar([EV.rate(1) REV.rate(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('PaleGreen'),'EdgeColor','none');hold on;
bar([0 0 EV.rate(2) REV.rate(2) 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 EV.rate(3) REV.rate(3) 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 EV.rate(4) REV.rate(4) 0 0 0 0 0 0 0 0],'FaceColor',rgb('SeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 EV.rate(5) REV.rate(5) 0 0 0 0 0 0],'FaceColor',rgb('LimeGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 EV.rate(6) REV.rate(6) 0 0 0 0],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(7) REV.rate(7) 0 0],'FaceColor',rgb('Green'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(8) REV.rate(8)],'FaceColor',rgb('DarkGreen'),'EdgeColor','none');hold on;
xlabel('EV/REVs per BLA cell FR octiles')
ylim([0 0.07])
subplot(3,2,2)
[h,bins]=hist(log(PyrBLAFR),100);
bar(bins,h,'EdgeColor','none','FaceColor',rgb('LightSteelBlue'));
hold on;
PlotHVLines(log(BLAFRoct),'v','r');

for i=1:8
  R_runpost=corrcoef(Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,3),Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,2));r_runpost.rate(i)=R_runpost(1,2);
  R_runpre=corrcoef(Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,1),Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,2));r_runpre.rate(i)=R_runpre(1,2);
  R_prepost=corrcoef(Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,1),Allcorr.sws(bothpyr&HpcFRgroups==i&nonan,3));r_prepost.rate(i)=R_prepost(1,2);
  EV.rate(i)=((r_runpost.rate(i)-r_runpre.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpre.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
  REV.rate(i)=((r_runpre.rate(i)-r_runpost.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpost.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
end  
 
 
subplot(3,2,3)
bar([EV.rate(1) REV.rate(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('PaleGreen'),'EdgeColor','none');hold on;
bar([0 0 EV.rate(2) REV.rate(2) 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 EV.rate(3) REV.rate(3) 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 EV.rate(4) REV.rate(4) 0 0 0 0 0 0 0 0],'FaceColor',rgb('SeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 EV.rate(5) REV.rate(5) 0 0 0 0 0 0],'FaceColor',rgb('LimeGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 EV.rate(6) REV.rate(6) 0 0 0 0],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(7) REV.rate(7) 0 0],'FaceColor',rgb('Green'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(8) REV.rate(8)],'FaceColor',rgb('DarkGreen'),'EdgeColor','none');hold on;
xlabel('EV/REVs per Hpc cell FR octiles')
ylim([0 0.07])
subplot(3,2,4)
[h,bins]=hist(log(PyrHpcFR),100);
bar(bins,h,'EdgeColor','none','FaceColor',rgb('LightSteelBlue'))
hold on;
PlotHVLines(log(HpcFRoct),'v','r')


for i=1:8
  R_runpost=corrcoef(Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,3),Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,2));r_runpost.rate(i)=R_runpost(1,2);
  R_runpre=corrcoef(Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,1),Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,2));r_runpre.rate(i)=R_runpre(1,2);
  R_prepost=corrcoef(Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,1),Allcorr.sws(bothpyr&MeanFRgroup==i&nonan,3));r_prepost.rate(i)=R_prepost(1,2);
  EV.rate(i)=((r_runpost.rate(i)-r_runpre.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpre.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
  REV.rate(i)=((r_runpre.rate(i)-r_runpost.rate(i)*r_prepost.rate(i))/sqrt((1-r_runpost.rate(i)^2)*(1-r_prepost.rate(i)^2)))^2;
end  

subplot(3,2,5);
bar([EV.rate(1) REV.rate(1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('PaleGreen'),'EdgeColor','none');hold on;
bar([0 0 EV.rate(2) REV.rate(2) 0 0 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('DarkSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 EV.rate(3) REV.rate(3) 0 0 0 0 0 0 0 0 0 0],'FaceColor',rgb('MediumSeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 EV.rate(4) REV.rate(4) 0 0 0 0 0 0 0 0],'FaceColor',rgb('SeaGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 EV.rate(5) REV.rate(5) 0 0 0 0 0 0],'FaceColor',rgb('LimeGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 EV.rate(6) REV.rate(6) 0 0 0 0],'FaceColor',rgb('ForestGreen'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(7) REV.rate(7) 0 0],'FaceColor',rgb('Green'),'EdgeColor','none');hold on;
bar([0 0 0 0 0 0 0 0 0 0 0 0 0 0 EV.rate(8) REV.rate(8)],'FaceColor',rgb('DarkGreen'),'EdgeColor','none');hold on;
xlabel('EV/REVs per Mean BLA/Hpc pair FR octiles');
ylim([0 0.07])
subplot(3,2,6);
[h,bins]=hist(log(PyrMeanFR),100);
bar(bins,h,'EdgeColor','none','FaceColor',rgb('LightSteelBlue'))
hold on;
PlotHVLines(log(MeanFRoct),'v','r')

%  figure('Position',[367 25 1014 964]);
%  subplot(3,3,1)
%  plot(Allcorr.run(ripup.runup,1),Allcorr.run(ripup.runup,3),'.','MarkerEdgeColor',rgb('Crimson'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  ylabel('Postrun corr. Ripple Upmod');
%  subplot(3,3,2)
%  plot(Allcorr.run(ripup.runnot,1),Allcorr.run(ripup.runnot,3),'.','MarkerEdgeColor',rgb('DarkGrey'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,3)
%  plot(Allcorr.run(ripup.rundown,1),Allcorr.run(ripup.rundown,3),'.','MarkerEdgeColor',rgb('FireBrick'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,4)
%  plot(Allcorr.run(ripnone.runup,1),Allcorr.run(ripnone.runup,3),'.','MarkerEdgeColor',rgb('SeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  ylabel('Postrun corr. Ripple Nonmod');
%  subplot(3,3,5)
%  plot(Allcorr.run(ripnone.runnot,1),Allcorr.run(ripnone.runnot,3),'.','MarkerEdgeColor',rgb('DarkSeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,6)
%  plot(Allcorr.run(ripnone.rundown,1),Allcorr.run(ripnone.rundown,3),'.','MarkerEdgeColor',rgb('MediumSeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,7)
%  plot(Allcorr.run(ripdown.runup,1),Allcorr.run(ripdown.runup,3),'.','MarkerEdgeColor',rgb('SteelBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN Sig. Pos. Corr.');
%  ylabel('Postrun corr. Ripple Downmod');
%  subplot(3,3,8)
%  plot(Allcorr.run(ripdown.runnot,1),Allcorr.run(ripdown.runnot,3),'.','MarkerEdgeColor',rgb('LightBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN Non-sig corr.')
%  subplot(3,3,9)
%  plot(Allcorr.run(ripdown.rundown,1),Allcorr.run(ripdown.rundown,3),'.','MarkerEdgeColor',rgb('DodgerBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN sig neg. corr. ')
%  
%  
%  figure('Position',[367 25 1014 964]);
%  subplot(3,3,1)
%  plot(Allcorr.run(ripup.runup,2),Allcorr.run(ripup.runup,3),'.','MarkerEdgeColor',rgb('Crimson'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  ylabel('Postrun corr. Ripple Upmod');
%  subplot(3,3,2)
%  plot(Allcorr.run(ripup.runnot,2),Allcorr.run(ripup.runnot,3),'.','MarkerEdgeColor',rgb('DarkGrey'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,3)
%  plot(Allcorr.run(ripup.rundown,2),Allcorr.run(ripup.rundown,3),'.','MarkerEdgeColor',rgb('FireBrick'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,4)
%  plot(Allcorr.run(ripnone.runup,2),Allcorr.run(ripnone.runup,3),'.','MarkerEdgeColor',rgb('SeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  ylabel('Postrun corr. Ripple Nonmod');
%  subplot(3,3,5)
%  plot(Allcorr.run(ripnone.runnot,2),Allcorr.run(ripnone.runnot,3),'.','MarkerEdgeColor',rgb('DarkSeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,6)
%  plot(Allcorr.run(ripnone.rundown,2),Allcorr.run(ripnone.rundown,3),'.','MarkerEdgeColor',rgb('MediumSeaGreen'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  subplot(3,3,7)
%  plot(Allcorr.run(ripdown.runup,2),Allcorr.run(ripdown.runup,3),'.','MarkerEdgeColor',rgb('SteelBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN Sig. Pos. Corr.');
%  ylabel('Postrun corr. Ripple Downmod');
%  subplot(3,3,8)
%  plot(Allcorr.run(ripdown.runnot,2),Allcorr.run(ripdown.runnot,3),'.','MarkerEdgeColor',rgb('LightBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN Non-sig corr.')
%  subplot(3,3,9)
%  plot(Allcorr.run(ripdown.rundown,2),Allcorr.run(ripdown.rundown,3),'.','MarkerEdgeColor',rgb('DodgerBlue'),'MarkerSize',15);
%  lsline;
%  xlim([-0.25 0.25]);ylim([-0.25 0.25]);
%  xlabel('PreRun corr. RUN sig neg. corr. ')
