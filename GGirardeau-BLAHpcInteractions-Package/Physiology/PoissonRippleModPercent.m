function [bla,hpc] = PoissonRippleModPercent(pval)
% PoissonRippleModPercent - Calculates percent of ripple-modulated cells in various structure
%
%  USAGE
%
%    [bla,hpc] = PoissonRippleModPercent (pval)
%
%    pval       pvalue under which a cell is considered modulated.
%
%  OUTPUT
%
%    bla    percentages for BLA
%    hpc    percentages for hippocampus
%
%  NOTES
%
%  SEE
%
%    See also : Eran's PoissonTest, PoissonRippleMod
%
% Jan 2014 by Gabrielle Girardeau calling Eran Stark's PoissonTest.m
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


load('/media/Data-01/All-Rats/PoissonRippleMod.mat');
load('/media/Data-01/All-Rats/AllRats-FinalType.mat');
load('/media/Data-01/All-Rats/SpikeParameters.mat');
load('/media/Data-01/All-Rats/Structures/structures.mat');

%  pval=0.001;
finaltype=finalType(ismember(finalType(:,1:4),poissonripplemod(:,1:4),'rows'),5);
ccgtype=SpikeParameters(ismember(SpikeParameters(:,1:4),poissonripplemod(:,1:4),'rows'),end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amygdala (BLA)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
bla.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),BLA,'rows'),:);
bla.finaltype=finaltype(ismember(poissonripplemod(:,1:3),BLA,'rows'));
bla.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),BLA,'rows'));

%  pval=0.05/size(bla.stats,1);

% counts
bla.nall.final=size(bla.stats,1);
bla.nall.ccg=sum(bla.ccgtype~=0);
bla.nint.final=sum(bla.finaltype==2);
bla.nint.ccg=sum(bla.ccgtype==2);
bla.npyr.final=sum(bla.finaltype==1);
bla.npyr.ccg=sum(bla.ccgtype==1);

bla.allsig=sum(bla.stats(:,6)<pval/2|bla.stats(:,7)<pval/2);
bla.allsigup=sum(bla.stats(:,6)<pval/2);
bla.allsigdown=sum(bla.stats(:,7)<pval/2);

% int final
bla.intsig.final=sum((bla.stats(:,6)<pval/2|bla.stats(:,7)<pval/2)&bla.finaltype==2);
bla.intsigup.final=sum(bla.stats(:,6)<pval/2&bla.finaltype==2);
bla.intsigdown.final=sum(bla.stats(:,7)<pval/2&bla.finaltype==2);
%int ccg
bla.intsig.ccg=sum((bla.stats(:,6)<pval/2|bla.stats(:,7)<pval/2)&bla.ccgtype==2);
bla.intsigup.ccg=sum(bla.stats(:,6)<pval/2&bla.ccgtype==2);
bla.intsigdown.ccg=sum(bla.stats(:,7)<pval/2&bla.ccgtype==2);

% pyr final
bla.pyrsig.final=sum((bla.stats(:,6)<pval/2|bla.stats(:,7)<pval/2)&bla.finaltype==1);
bla.pyrsigup.final=sum(bla.stats(:,6)<pval/2&bla.finaltype==1);
bla.pyrsigdown.final=sum(bla.stats(:,7)<pval/2&bla.finaltype==1);
% pyr ccg
bla.pyrsig.ccg=sum((bla.stats(:,6)<pval/2|bla.stats(:,7)<pval/2)&bla.ccgtype==1);
bla.pyrsigup.ccg=sum(bla.stats(:,6)<pval/2&bla.ccgtype==1);
bla.pyrsigdown.ccg=sum(bla.stats(:,7)<pval/2&bla.ccgtype==1);

% percents
bla.percentall=bla.allsig/bla.nall.final*100;
bla.percentallup=bla.allsigup/bla.nall.final*100;
bla.percentalldown=bla.allsigdown/bla.nall.final*100;

bla.percentint.final=bla.intsig.final/bla.nint.final*100;
bla.percentintup.final=bla.intsigup.final/bla.nint.final*100;
bla.percentintdown.final=bla.intsigdown.final/bla.nint.final*100;

bla.percentint.ccg=bla.intsig.ccg/bla.nint.ccg*100;
bla.percentintup.ccg=bla.intsigup.ccg/bla.nint.ccg*100;
bla.percentintdown.ccg=bla.intsigdown.ccg/bla.nint.ccg*100;

bla.percentpyr.final=bla.pyrsig.final/bla.npyr.final*100;
bla.percentpyrup.final=bla.pyrsigup.final/bla.npyr.final*100;
bla.percentpyrdown.final=bla.pyrsigdown.final/bla.npyr.final*100;

bla.percentpyr.ccg=bla.pyrsig.ccg/bla.npyr.ccg*100;
bla.percentpyrup.ccg=bla.pyrsigup.ccg/bla.npyr.ccg*100;
bla.percentpyrdown.ccg=bla.pyrsigdown.ccg/bla.npyr.ccg*100;

figure;
subplot(2,6,1);
bar([bla.percentalldown*(-1) 0 bla.percentpyrdown.final*(-1) bla.percentintdown.final*(-1) 0 bla.percentpyrdown.ccg*(-1) bla.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([bla.percentallup 0 bla.percentpyrup.final bla.percentintup.final 0 bla.percentpyrup.ccg bla.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'BLA, Poisson Ripple Mod';['nAll=' int2str(bla.nall.final) ' nPyr.f=' int2str(bla.npyr.final) ' nInt.f=' int2str(bla.nint.final)];[' nPyr.ccg=' int2str(bla.npyr.ccg) ' nInt.ccg=' int2str(bla.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amygdala (ALL - basomediolat)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
amy.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),basal,'rows'),:);
amy.finaltype=finaltype(ismember(poissonripplemod(:,1:3),basal,'rows'));
amy.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),basal,'rows'));

% counts
amy.nall.final=size(amy.stats,1);
amy.nall.ccg=sum(amy.ccgtype~=0);
amy.nint.final=sum(amy.finaltype==2);
amy.nint.ccg=sum(amy.ccgtype==2);
amy.npyr.final=sum(amy.finaltype==1);
amy.npyr.ccg=sum(amy.ccgtype==1);

amy.allsig=sum(amy.stats(:,6)<pval/2|amy.stats(:,7)<pval/2);
amy.allsigup=sum(amy.stats(:,6)<pval/2);
amy.allsigdown=sum(amy.stats(:,7)<pval/2);

% int final
amy.intsig.final=sum((amy.stats(:,6)<pval/2|amy.stats(:,7)<pval/2)&amy.finaltype==2);
amy.intsigup.final=sum(amy.stats(:,6)<pval/2&amy.finaltype==2);
amy.intsigdown.final=sum(amy.stats(:,7)<pval/2&amy.finaltype==2);
%int ccg
amy.intsig.ccg=sum((amy.stats(:,6)<pval/2|amy.stats(:,7)<pval/2)&amy.ccgtype==2);
amy.intsigup.ccg=sum(amy.stats(:,6)<pval/2&amy.ccgtype==2);
amy.intsigdown.ccg=sum(amy.stats(:,7)<pval/2&amy.ccgtype==2);

% pyr final
amy.pyrsig.final=sum((amy.stats(:,6)<pval/2|amy.stats(:,7)<pval/2)&amy.finaltype==1);
amy.pyrsigup.final=sum(amy.stats(:,6)<pval/2&amy.finaltype==1);
amy.pyrsigdown.final=sum(amy.stats(:,7)<pval/2&amy.finaltype==1);
% pyr ccg
amy.pyrsig.ccg=sum((amy.stats(:,6)<pval/2|amy.stats(:,7)<pval/2)&amy.ccgtype==1);
amy.pyrsigup.ccg=sum(amy.stats(:,6)<pval/2&amy.ccgtype==1);
amy.pyrsigdown.ccg=sum(amy.stats(:,7)<pval/2&amy.ccgtype==1);

% percents
amy.percentall=amy.allsig/amy.nall.final*100;
amy.percentallup=amy.allsigup/amy.nall.final*100;
amy.percentalldown=amy.allsigdown/amy.nall.final*100;

amy.percentint.final=amy.intsig.final/amy.nint.final*100;
amy.percentintup.final=amy.intsigup.final/amy.nint.final*100;
amy.percentintdown.final=amy.intsigdown.final/amy.nint.final*100;

amy.percentint.ccg=amy.intsig.ccg/amy.nint.ccg*100;
amy.percentintup.ccg=amy.intsigup.ccg/amy.nint.ccg*100;
amy.percentintdown.ccg=amy.intsigdown.ccg/amy.nint.ccg*100;

amy.percentpyr.final=amy.pyrsig.final/amy.npyr.final*100;
amy.percentpyrup.final=amy.pyrsigup.final/amy.npyr.final*100;
amy.percentpyrdown.final=amy.pyrsigdown.final/amy.npyr.final*100;

amy.percentpyr.ccg=amy.pyrsig.ccg/amy.npyr.ccg*100;
amy.percentpyrup.ccg=amy.pyrsigup.ccg/amy.npyr.ccg*100;
amy.percentpyrdown.ccg=amy.pyrsigdown.ccg/amy.npyr.ccg*100;

%  figure;
subplot(2,6,2);
bar([amy.percentalldown*(-1) 0 amy.percentpyrdown.final*(-1) amy.percentintdown.final*(-1) 0 amy.percentpyrdown.ccg*(-1) amy.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([amy.percentallup 0 amy.percentpyrup.final amy.percentintup.final 0 amy.percentpyrup.ccg amy.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'All B/L/M amy, Poisson Ripple Mod';['nAll=' int2str(amy.nall.final) ' nPyr.f=' int2str(amy.npyr.final) ' nInt.f=' int2str(amy.nint.final)];[' nPyr.ccg=' int2str(amy.npyr.ccg) ' nInt.ccg=' int2str(amy.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hippocampus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
hpc.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),Hpc,'rows'),:);
hpc.finaltype=finaltype(ismember(poissonripplemod(:,1:3),Hpc,'rows'));
hpc.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),Hpc,'rows'));

% counts
hpc.nall.final=size(hpc.stats,1);
hpc.nall.ccg=sum(hpc.ccgtype~=0);
hpc.nint.final=sum(hpc.finaltype==2);
hpc.nint.ccg=sum(hpc.ccgtype==2);
hpc.npyr.final=sum(hpc.finaltype==1);
hpc.npyr.ccg=sum(hpc.ccgtype==1);

hpc.allsig=sum(hpc.stats(:,6)<pval/2|hpc.stats(:,7)<pval/2);
hpc.allsigup=sum(hpc.stats(:,6)<pval/2);
hpc.allsigdown=sum(hpc.stats(:,7)<pval/2);

% int final
hpc.intsig.final=sum((hpc.stats(:,6)<pval/2|hpc.stats(:,7)<pval/2)&hpc.finaltype==2);
hpc.intsigup.final=sum(hpc.stats(:,6)<pval/2&hpc.finaltype==2);
hpc.intsigdown.final=sum(hpc.stats(:,7)<pval/2&hpc.finaltype==2);
%int ccg
hpc.intsig.ccg=sum((hpc.stats(:,6)<pval/2|hpc.stats(:,7)<pval/2)&hpc.ccgtype==2);
hpc.intsigup.ccg=sum(hpc.stats(:,6)<pval/2&hpc.ccgtype==2);
hpc.intsigdown.ccg=sum(hpc.stats(:,7)<pval/2&hpc.ccgtype==2);

% pyr final
hpc.pyrsig.final=sum((hpc.stats(:,6)<pval/2|hpc.stats(:,7)<pval/2)&hpc.finaltype==1);
hpc.pyrsigup.final=sum(hpc.stats(:,6)<pval/2&hpc.finaltype==1);
hpc.pyrsigdown.final=sum(hpc.stats(:,7)<pval/2&hpc.finaltype==1);
% pyr ccg
hpc.pyrsig.ccg=sum((hpc.stats(:,6)<pval/2|hpc.stats(:,7)<pval/2)&hpc.ccgtype==1);
hpc.pyrsigup.ccg=sum(hpc.stats(:,6)<pval/2&hpc.ccgtype==1);
hpc.pyrsigdown.ccg=sum(hpc.stats(:,7)<pval/2&hpc.ccgtype==1);

% percents
hpc.percentall=hpc.allsig/hpc.nall.final*100;
hpc.percentallup=hpc.allsigup/hpc.nall.final*100;
hpc.percentalldown=hpc.allsigdown/hpc.nall.final*100;

hpc.percentint.final=hpc.intsig.final/hpc.nint.final*100;
hpc.percentintup.final=hpc.intsigup.final/hpc.nint.final*100;
hpc.percentintdown.final=hpc.intsigdown.final/hpc.nint.final*100;

hpc.percentint.ccg=hpc.intsig.ccg/hpc.nint.ccg*100;
hpc.percentintup.ccg=hpc.intsigup.ccg/hpc.nint.ccg*100;
hpc.percentintdown.ccg=hpc.intsigdown.ccg/hpc.nint.ccg*100;

hpc.percentpyr.final=hpc.pyrsig.final/hpc.npyr.final*100;
hpc.percentpyrup.final=hpc.pyrsigup.final/hpc.npyr.final*100;
hpc.percentpyrdown.final=hpc.pyrsigdown.final/hpc.npyr.final*100;

hpc.percentpyr.ccg=hpc.pyrsig.ccg/hpc.npyr.ccg*100;
hpc.percentpyrup.ccg=hpc.pyrsigup.ccg/hpc.npyr.ccg*100;
hpc.percentpyrdown.ccg=hpc.pyrsigdown.ccg/hpc.npyr.ccg*100;

%  figure;
subplot(2,6,3);
bar([hpc.percentalldown*(-1) 0 hpc.percentpyrdown.final*(-1) hpc.percentintdown.final*(-1) 0 hpc.percentpyrdown.ccg*(-1) hpc.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([hpc.percentallup 0 hpc.percentpyrup.final hpc.percentintup.final 0 hpc.percentpyrup.ccg hpc.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'Hpc, Poisson Ripple Mod';['nAll=' int2str(hpc.nall.final) ' nPyr.f=' int2str(hpc.npyr.final) ' nInt.f=' int2str(hpc.nint.final)];[' nPyr.ccg=' int2str(hpc.npyr.ccg) ' nInt.ccg=' int2str(hpc.nint.ccg) ' p<' num2str(pval)]});
ylim([-10 100]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
pir.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),Pir,'rows'),:);
pir.finaltype=finaltype(ismember(poissonripplemod(:,1:3),Pir,'rows'));
pir.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),Pir,'rows'));

% counts
pir.nall.final=size(pir.stats,1);
pir.nall.ccg=sum(pir.ccgtype~=0);
pir.nint.final=sum(pir.finaltype==2);
pir.nint.ccg=sum(pir.ccgtype==2);
pir.npyr.final=sum(pir.finaltype==1);
pir.npyr.ccg=sum(pir.ccgtype==1);

pir.allsig=sum(pir.stats(:,6)<pval/2|pir.stats(:,7)<pval/2);
pir.allsigup=sum(pir.stats(:,6)<pval/2);
pir.allsigdown=sum(pir.stats(:,7)<pval/2);

% int final
pir.intsig.final=sum((pir.stats(:,6)<pval/2|pir.stats(:,7)<pval/2)&pir.finaltype==2);
pir.intsigup.final=sum(pir.stats(:,6)<pval/2&pir.finaltype==2);
pir.intsigdown.final=sum(pir.stats(:,7)<pval/2&pir.finaltype==2);
%int ccg
pir.intsig.ccg=sum((pir.stats(:,6)<pval/2|pir.stats(:,7)<pval/2)&pir.ccgtype==2);
pir.intsigup.ccg=sum(pir.stats(:,6)<pval/2&pir.ccgtype==2);
pir.intsigdown.ccg=sum(pir.stats(:,7)<pval/2&pir.ccgtype==2);

% pyr final
pir.pyrsig.final=sum((pir.stats(:,6)<pval/2|pir.stats(:,7)<pval/2)&pir.finaltype==1);
pir.pyrsigup.final=sum(pir.stats(:,6)<pval/2&pir.finaltype==1);
pir.pyrsigdown.final=sum(pir.stats(:,7)<pval/2&pir.finaltype==1);
% pyr ccg
pir.pyrsig.ccg=sum((pir.stats(:,6)<pval/2|pir.stats(:,7)<pval/2)&pir.ccgtype==1);
pir.pyrsigup.ccg=sum(pir.stats(:,6)<pval/2&pir.ccgtype==1);
pir.pyrsigdown.ccg=sum(pir.stats(:,7)<pval/2&pir.ccgtype==1);

% percents
pir.percentall=pir.allsig/pir.nall.final*100;
pir.percentallup=pir.allsigup/pir.nall.final*100;
pir.percentalldown=pir.allsigdown/pir.nall.final*100;

pir.percentint.final=pir.intsig.final/pir.nint.final*100;
pir.percentintup.final=pir.intsigup.final/pir.nint.final*100;
pir.percentintdown.final=pir.intsigdown.final/pir.nint.final*100;

pir.percentint.ccg=pir.intsig.ccg/pir.nint.ccg*100;
pir.percentintup.ccg=pir.intsigup.ccg/pir.nint.ccg*100;
pir.percentintdown.ccg=pir.intsigdown.ccg/pir.nint.ccg*100;

pir.percentpyr.final=pir.pyrsig.final/pir.npyr.final*100;
pir.percentpyrup.final=pir.pyrsigup.final/pir.npyr.final*100;
pir.percentpyrdown.final=pir.pyrsigdown.final/pir.npyr.final*100;

pir.percentpyr.ccg=pir.pyrsig.ccg/pir.npyr.ccg*100;
pir.percentpyrup.ccg=pir.pyrsigup.ccg/pir.npyr.ccg*100;
pir.percentpyrdown.ccg=pir.pyrsigdown.ccg/pir.npyr.ccg*100;

%  figure;
subplot(2,6,4);
bar([pir.percentalldown*(-1) 0 pir.percentpyrdown.final*(-1) pir.percentintdown.final*(-1) 0 pir.percentpyrdown.ccg*(-1) pir.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([pir.percentallup 0 pir.percentpyrup.final pir.percentintup.final 0 pir.percentpyrup.ccg pir.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'Pir, Poisson Ripple Mod';['nAll=' int2str(pir.nall.final) ' nPyr.f=' int2str(pir.npyr.final) ' nInt.f=' int2str(pir.nint.final)];[' nPyr.ccg=' int2str(pir.npyr.ccg) ' nInt.ccg=' int2str(pir.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Olfact (Pir + DEn +VEn)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
olf.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),olfact,'rows'),:);
olf.finaltype=finaltype(ismember(poissonripplemod(:,1:3),olfact,'rows'));
olf.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),olfact,'rows'));

% counts
olf.nall.final=size(olf.stats,1);
olf.nall.ccg=sum(olf.ccgtype~=0);
olf.nint.final=sum(olf.finaltype==2);
olf.nint.ccg=sum(olf.ccgtype==2);
olf.npyr.final=sum(olf.finaltype==1);
olf.npyr.ccg=sum(olf.ccgtype==1);

olf.allsig=sum(olf.stats(:,6)<pval/2|olf.stats(:,7)<pval/2);
olf.allsigup=sum(olf.stats(:,6)<pval/2);
olf.allsigdown=sum(olf.stats(:,7)<pval/2);

% int final
olf.intsig.final=sum((olf.stats(:,6)<pval/2|olf.stats(:,7)<pval/2)&olf.finaltype==2);
olf.intsigup.final=sum(olf.stats(:,6)<pval/2&olf.finaltype==2);
olf.intsigdown.final=sum(olf.stats(:,7)<pval/2&olf.finaltype==2);
%int ccg
olf.intsig.ccg=sum((olf.stats(:,6)<pval/2|olf.stats(:,7)<pval/2)&olf.ccgtype==2);
olf.intsigup.ccg=sum(olf.stats(:,6)<pval/2&olf.ccgtype==2);
olf.intsigdown.ccg=sum(olf.stats(:,7)<pval/2&olf.ccgtype==2);

% pyr final
olf.pyrsig.final=sum((olf.stats(:,6)<pval/2|olf.stats(:,7)<pval/2)&olf.finaltype==1);
olf.pyrsigup.final=sum(olf.stats(:,6)<pval/2&olf.finaltype==1);
olf.pyrsigdown.final=sum(olf.stats(:,7)<pval/2&olf.finaltype==1);
% pyr ccg
olf.pyrsig.ccg=sum((olf.stats(:,6)<pval/2|olf.stats(:,7)<pval/2)&olf.ccgtype==1);
olf.pyrsigup.ccg=sum(olf.stats(:,6)<pval/2&olf.ccgtype==1);
olf.pyrsigdown.ccg=sum(olf.stats(:,7)<pval/2&olf.ccgtype==1);

% percents
olf.percentall=olf.allsig/olf.nall.final*100;
olf.percentallup=olf.allsigup/olf.nall.final*100;
olf.percentalldown=olf.allsigdown/olf.nall.final*100;

olf.percentint.final=olf.intsig.final/olf.nint.final*100;
olf.percentintup.final=olf.intsigup.final/olf.nint.final*100;
olf.percentintdown.final=olf.intsigdown.final/olf.nint.final*100;

olf.percentint.ccg=olf.intsig.ccg/olf.nint.ccg*100;
olf.percentintup.ccg=olf.intsigup.ccg/olf.nint.ccg*100;
olf.percentintdown.ccg=olf.intsigdown.ccg/olf.nint.ccg*100;

olf.percentpyr.final=olf.pyrsig.final/olf.npyr.final*100;
olf.percentpyrup.final=olf.pyrsigup.final/olf.npyr.final*100;
olf.percentpyrdown.final=olf.pyrsigdown.final/olf.npyr.final*100;

olf.percentpyr.ccg=olf.pyrsig.ccg/olf.npyr.ccg*100;
olf.percentpyrup.ccg=olf.pyrsigup.ccg/olf.npyr.ccg*100;
olf.percentpyrdown.ccg=olf.pyrsigdown.ccg/olf.npyr.ccg*100;

%  figure;
subplot(2,6,5);
bar([olf.percentalldown*(-1) 0 olf.percentpyrdown.final*(-1) olf.percentintdown.final*(-1) 0 olf.percentpyrdown.ccg*(-1) olf.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([olf.percentallup 0 olf.percentpyrup.final olf.percentintup.final 0 olf.percentpyrup.ccg olf.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'Pir + D/VEn, Poisson Ripple Mod';['nAll=' int2str(olf.nall.final) ' nPyr.f=' int2str(olf.npyr.final) ' nInt.f=' int2str(olf.nint.final)];[' nPyr.ccg=' int2str(olf.npyr.ccg) ' nInt.ccg=' int2str(olf.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CeCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
cecm.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),CeCM,'rows'),:);
cecm.finaltype=finaltype(ismember(poissonripplemod(:,1:3),CeCM,'rows'));
cecm.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),CeCM,'rows'));

% counts
cecm.nall.final=size(cecm.stats,1);
cecm.nall.ccg=sum(cecm.ccgtype~=0);
cecm.nint.final=sum(cecm.finaltype==2);
cecm.nint.ccg=sum(cecm.ccgtype==2);
cecm.npyr.final=sum(cecm.finaltype==1);
cecm.npyr.ccg=sum(cecm.ccgtype==1);

cecm.allsig=sum(cecm.stats(:,6)<pval/2|cecm.stats(:,7)<pval/2);
cecm.allsigup=sum(cecm.stats(:,6)<pval/2);
cecm.allsigdown=sum(cecm.stats(:,7)<pval/2);

% int final
cecm.intsig.final=sum((cecm.stats(:,6)<pval/2|cecm.stats(:,7)<pval/2)&cecm.finaltype==2);
cecm.intsigup.final=sum(cecm.stats(:,6)<pval/2&cecm.finaltype==2);
cecm.intsigdown.final=sum(cecm.stats(:,7)<pval/2&cecm.finaltype==2);
%int ccg
cecm.intsig.ccg=sum((cecm.stats(:,6)<pval/2|cecm.stats(:,7)<pval/2)&cecm.ccgtype==2);
cecm.intsigup.ccg=sum(cecm.stats(:,6)<pval/2&cecm.ccgtype==2);
cecm.intsigdown.ccg=sum(cecm.stats(:,7)<pval/2&cecm.ccgtype==2);

% pyr final
cecm.pyrsig.final=sum((cecm.stats(:,6)<pval/2|cecm.stats(:,7)<pval/2)&cecm.finaltype==1);
cecm.pyrsigup.final=sum(cecm.stats(:,6)<pval/2&cecm.finaltype==1);
cecm.pyrsigdown.final=sum(cecm.stats(:,7)<pval/2&cecm.finaltype==1);
% pyr ccg
cecm.pyrsig.ccg=sum((cecm.stats(:,6)<pval/2|cecm.stats(:,7)<pval/2)&cecm.ccgtype==1);
cecm.pyrsigup.ccg=sum(cecm.stats(:,6)<pval/2&cecm.ccgtype==1);
cecm.pyrsigdown.ccg=sum(cecm.stats(:,7)<pval/2&cecm.ccgtype==1);

% percents
cecm.percentall=cecm.allsig/cecm.nall.final*100;
cecm.percentallup=cecm.allsigup/cecm.nall.final*100;
cecm.percentalldown=cecm.allsigdown/cecm.nall.final*100;

cecm.percentint.final=cecm.intsig.final/cecm.nint.final*100;
cecm.percentintup.final=cecm.intsigup.final/cecm.nint.final*100;
cecm.percentintdown.final=cecm.intsigdown.final/cecm.nint.final*100;

cecm.percentint.ccg=cecm.intsig.ccg/cecm.nint.ccg*100;
cecm.percentintup.ccg=cecm.intsigup.ccg/cecm.nint.ccg*100;
cecm.percentintdown.ccg=cecm.intsigdown.ccg/cecm.nint.ccg*100;

cecm.percentpyr.final=cecm.pyrsig.final/cecm.npyr.final*100;
cecm.percentpyrup.final=cecm.pyrsigup.final/cecm.npyr.final*100;
cecm.percentpyrdown.final=cecm.pyrsigdown.final/cecm.npyr.final*100;

cecm.percentpyr.ccg=cecm.pyrsig.ccg/cecm.npyr.ccg*100;
cecm.percentpyrup.ccg=cecm.pyrsigup.ccg/cecm.npyr.ccg*100;
cecm.percentpyrdown.ccg=cecm.pyrsigdown.ccg/cecm.npyr.ccg*100;

%  figure;
subplot(2,6,6);
bar([cecm.percentalldown*(-1) 0 cecm.percentpyrdown.final*(-1) cecm.percentintdown.final*(-1) 0 cecm.percentpyrdown.ccg*(-1) cecm.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([cecm.percentallup 0 cecm.percentpyrup.final cecm.percentintup.final 0 cecm.percentpyrup.ccg cecm.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'CeCM, Poisson Ripple Mod';['nAll=' int2str(cecm.nall.final) ' nPyr.f=' int2str(cecm.npyr.final) ' nInt.f=' int2str(cecm.nint.final)];[' nPyr.ccg=' int2str(cecm.npyr.ccg) ' nInt.ccg=' int2str(cecm.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LaDL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
ladl.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),LaDL,'rows'),:);
ladl.finaltype=finaltype(ismember(poissonripplemod(:,1:3),LaDL,'rows'));
ladl.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),LaDL,'rows'));

% counts
ladl.nall.final=size(ladl.stats,1);
ladl.nall.ccg=sum(ladl.ccgtype~=0);
ladl.nint.final=sum(ladl.finaltype==2);
ladl.nint.ccg=sum(ladl.ccgtype==2);
ladl.npyr.final=sum(ladl.finaltype==1);
ladl.npyr.ccg=sum(ladl.ccgtype==1);

ladl.allsig=sum(ladl.stats(:,6)<pval/2|ladl.stats(:,7)<pval/2);
ladl.allsigup=sum(ladl.stats(:,6)<pval/2);
ladl.allsigdown=sum(ladl.stats(:,7)<pval/2);

% int final
ladl.intsig.final=sum((ladl.stats(:,6)<pval/2|ladl.stats(:,7)<pval/2)&ladl.finaltype==2);
ladl.intsigup.final=sum(ladl.stats(:,6)<pval/2&ladl.finaltype==2);
ladl.intsigdown.final=sum(ladl.stats(:,7)<pval/2&ladl.finaltype==2);
%int ccg
ladl.intsig.ccg=sum((ladl.stats(:,6)<pval/2|ladl.stats(:,7)<pval/2)&ladl.ccgtype==2);
ladl.intsigup.ccg=sum(ladl.stats(:,6)<pval/2&ladl.ccgtype==2);
ladl.intsigdown.ccg=sum(ladl.stats(:,7)<pval/2&ladl.ccgtype==2);

% pyr final
ladl.pyrsig.final=sum((ladl.stats(:,6)<pval/2|ladl.stats(:,7)<pval/2)&ladl.finaltype==1);
ladl.pyrsigup.final=sum(ladl.stats(:,6)<pval/2&ladl.finaltype==1);
ladl.pyrsigdown.final=sum(ladl.stats(:,7)<pval/2&ladl.finaltype==1);
% pyr ccg
ladl.pyrsig.ccg=sum((ladl.stats(:,6)<pval/2|ladl.stats(:,7)<pval/2)&ladl.ccgtype==1);
ladl.pyrsigup.ccg=sum(ladl.stats(:,6)<pval/2&ladl.ccgtype==1);
ladl.pyrsigdown.ccg=sum(ladl.stats(:,7)<pval/2&ladl.ccgtype==1);

% percents
ladl.percentall=ladl.allsig/ladl.nall.final*100;
ladl.percentallup=ladl.allsigup/ladl.nall.final*100;
ladl.percentalldown=ladl.allsigdown/ladl.nall.final*100;

ladl.percentint.final=ladl.intsig.final/ladl.nint.final*100;
ladl.percentintup.final=ladl.intsigup.final/ladl.nint.final*100;
ladl.percentintdown.final=ladl.intsigdown.final/ladl.nint.final*100;

ladl.percentint.ccg=ladl.intsig.ccg/ladl.nint.ccg*100;
ladl.percentintup.ccg=ladl.intsigup.ccg/ladl.nint.ccg*100;
ladl.percentintdown.ccg=ladl.intsigdown.ccg/ladl.nint.ccg*100;

ladl.percentpyr.final=ladl.pyrsig.final/ladl.npyr.final*100;
ladl.percentpyrup.final=ladl.pyrsigup.final/ladl.npyr.final*100;
ladl.percentpyrdown.final=ladl.pyrsigdown.final/ladl.npyr.final*100;

ladl.percentpyr.ccg=ladl.pyrsig.ccg/ladl.npyr.ccg*100;
ladl.percentpyrup.ccg=ladl.pyrsigup.ccg/ladl.npyr.ccg*100;
ladl.percentpyrdown.ccg=ladl.pyrsigdown.ccg/ladl.npyr.ccg*100;

%  figure;
subplot(2,6,7);
bar([ladl.percentalldown*(-1) 0 ladl.percentpyrdown.final*(-1) ladl.percentintdown.final*(-1) 0 ladl.percentpyrdown.ccg*(-1) ladl.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([ladl.percentallup 0 ladl.percentpyrup.final ladl.percentintup.final 0 ladl.percentpyrup.ccg ladl.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'LaDL, Poisson Ripple Mod';['nAll=' int2str(ladl.nall.final) ' nPyr.f=' int2str(ladl.npyr.final) ' nInt.f=' int2str(ladl.nint.final)];[' nPyr.ccg=' int2str(ladl.npyr.ccg) ' nInt.ccg=' int2str(ladl.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BMA + BMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
bma.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),[BMA;BMP],'rows'),:);
bma.finaltype=finaltype(ismember(poissonripplemod(:,1:3),[BMA;BMP],'rows'));
bma.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),[BMA;BMP],'rows'));

% counts
bma.nall.final=size(bma.stats,1);
bma.nall.ccg=sum(bma.ccgtype~=0);
bma.nint.final=sum(bma.finaltype==2);
bma.nint.ccg=sum(bma.ccgtype==2);
bma.npyr.final=sum(bma.finaltype==1);
bma.npyr.ccg=sum(bma.ccgtype==1);

bma.allsig=sum(bma.stats(:,6)<pval/2|bma.stats(:,7)<pval/2);
bma.allsigup=sum(bma.stats(:,6)<pval/2);
bma.allsigdown=sum(bma.stats(:,7)<pval/2);

% int final
bma.intsig.final=sum((bma.stats(:,6)<pval/2|bma.stats(:,7)<pval/2)&bma.finaltype==2);
bma.intsigup.final=sum(bma.stats(:,6)<pval/2&bma.finaltype==2);
bma.intsigdown.final=sum(bma.stats(:,7)<pval/2&bma.finaltype==2);
%int ccg
bma.intsig.ccg=sum((bma.stats(:,6)<pval/2|bma.stats(:,7)<pval/2)&bma.ccgtype==2);
bma.intsigup.ccg=sum(bma.stats(:,6)<pval/2&bma.ccgtype==2);
bma.intsigdown.ccg=sum(bma.stats(:,7)<pval/2&bma.ccgtype==2);

% pyr final
bma.pyrsig.final=sum((bma.stats(:,6)<pval/2|bma.stats(:,7)<pval/2)&bma.finaltype==1);
bma.pyrsigup.final=sum(bma.stats(:,6)<pval/2&bma.finaltype==1);
bma.pyrsigdown.final=sum(bma.stats(:,7)<pval/2&bma.finaltype==1);
% pyr ccg
bma.pyrsig.ccg=sum((bma.stats(:,6)<pval/2|bma.stats(:,7)<pval/2)&bma.ccgtype==1);
bma.pyrsigup.ccg=sum(bma.stats(:,6)<pval/2&bma.ccgtype==1);
bma.pyrsigdown.ccg=sum(bma.stats(:,7)<pval/2&bma.ccgtype==1);

% percents
bma.percentall=bma.allsig/bma.nall.final*100;
bma.percentallup=bma.allsigup/bma.nall.final*100;
bma.percentalldown=bma.allsigdown/bma.nall.final*100;

bma.percentint.final=bma.intsig.final/bma.nint.final*100;
bma.percentintup.final=bma.intsigup.final/bma.nint.final*100;
bma.percentintdown.final=bma.intsigdown.final/bma.nint.final*100;

bma.percentint.ccg=bma.intsig.ccg/bma.nint.ccg*100;
bma.percentintup.ccg=bma.intsigup.ccg/bma.nint.ccg*100;
bma.percentintdown.ccg=bma.intsigdown.ccg/bma.nint.ccg*100;

bma.percentpyr.final=bma.pyrsig.final/bma.npyr.final*100;
bma.percentpyrup.final=bma.pyrsigup.final/bma.npyr.final*100;
bma.percentpyrdown.final=bma.pyrsigdown.final/bma.npyr.final*100;

bma.percentpyr.ccg=bma.pyrsig.ccg/bma.npyr.ccg*100;
bma.percentpyrup.ccg=bma.pyrsigup.ccg/bma.npyr.ccg*100;
bma.percentpyrdown.ccg=bma.pyrsigdown.ccg/bma.npyr.ccg*100;

%  figure;
subplot(2,6,8);
bar([bma.percentalldown*(-1) 0 bma.percentpyrdown.final*(-1) bma.percentintdown.final*(-1) 0 bma.percentpyrdown.ccg*(-1) bma.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([bma.percentallup 0 bma.percentpyrup.final bma.percentintup.final 0 bma.percentpyrup.ccg bma.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'BMA+BMP, Poisson Ripple Mod';['nAll=' int2str(bma.nall.final) ' nPyr.f=' int2str(bma.npyr.final) ' nInt.f=' int2str(bma.nint.final)];[' nPyr.ccg=' int2str(bma.npyr.ccg) ' nInt.ccg=' int2str(bma.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
blv.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),BLV,'rows'),:);
blv.finaltype=finaltype(ismember(poissonripplemod(:,1:3),BLV,'rows'));
blv.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),BLV,'rows'));

% counts
blv.nall.final=size(blv.stats,1);
blv.nall.ccg=sum(blv.ccgtype~=0);
blv.nint.final=sum(blv.finaltype==2);
blv.nint.ccg=sum(blv.ccgtype==2);
blv.npyr.final=sum(blv.finaltype==1);
blv.npyr.ccg=sum(blv.ccgtype==1);

blv.allsig=sum(blv.stats(:,6)<pval/2|blv.stats(:,7)<pval/2);
blv.allsigup=sum(blv.stats(:,6)<pval/2);
blv.allsigdown=sum(blv.stats(:,7)<pval/2);

% int final
blv.intsig.final=sum((blv.stats(:,6)<pval/2|blv.stats(:,7)<pval/2)&blv.finaltype==2);
blv.intsigup.final=sum(blv.stats(:,6)<pval/2&blv.finaltype==2);
blv.intsigdown.final=sum(blv.stats(:,7)<pval/2&blv.finaltype==2);
%int ccg
blv.intsig.ccg=sum((blv.stats(:,6)<pval/2|blv.stats(:,7)<pval/2)&blv.ccgtype==2);
blv.intsigup.ccg=sum(blv.stats(:,6)<pval/2&blv.ccgtype==2);
blv.intsigdown.ccg=sum(blv.stats(:,7)<pval/2&blv.ccgtype==2);

% pyr final
blv.pyrsig.final=sum((blv.stats(:,6)<pval/2|blv.stats(:,7)<pval/2)&blv.finaltype==1);
blv.pyrsigup.final=sum(blv.stats(:,6)<pval/2&blv.finaltype==1);
blv.pyrsigdown.final=sum(blv.stats(:,7)<pval/2&blv.finaltype==1);
% pyr ccg
blv.pyrsig.ccg=sum((blv.stats(:,6)<pval/2|blv.stats(:,7)<pval/2)&blv.ccgtype==1);
blv.pyrsigup.ccg=sum(blv.stats(:,6)<pval/2&blv.ccgtype==1);
blv.pyrsigdown.ccg=sum(blv.stats(:,7)<pval/2&blv.ccgtype==1);

% percents
blv.percentall=blv.allsig/blv.nall.final*100;
blv.percentallup=blv.allsigup/blv.nall.final*100;
blv.percentalldown=blv.allsigdown/blv.nall.final*100;

blv.percentint.final=blv.intsig.final/blv.nint.final*100;
blv.percentintup.final=blv.intsigup.final/blv.nint.final*100;
blv.percentintdown.final=blv.intsigdown.final/blv.nint.final*100;

blv.percentint.ccg=blv.intsig.ccg/blv.nint.ccg*100;
blv.percentintup.ccg=blv.intsigup.ccg/blv.nint.ccg*100;
blv.percentintdown.ccg=blv.intsigdown.ccg/blv.nint.ccg*100;

blv.percentpyr.final=blv.pyrsig.final/blv.npyr.final*100;
blv.percentpyrup.final=blv.pyrsigup.final/blv.npyr.final*100;
blv.percentpyrdown.final=blv.pyrsigdown.final/blv.npyr.final*100;

blv.percentpyr.ccg=blv.pyrsig.ccg/blv.npyr.ccg*100;
blv.percentpyrup.ccg=blv.pyrsigup.ccg/blv.npyr.ccg*100;
blv.percentpyrdown.ccg=blv.pyrsigdown.ccg/blv.npyr.ccg*100;

%  figure;
subplot(2,6,10);
bar([blv.percentalldown*(-1) 0 blv.percentpyrdown.final*(-1) blv.percentintdown.final*(-1) 0 blv.percentpyrdown.ccg*(-1) blv.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([blv.percentallup 0 blv.percentpyrup.final blv.percentintup.final 0 blv.percentpyrup.ccg blv.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'BLV, Poisson Ripple Mod';['nAll=' int2str(blv.nall.final) ' nPyr.f=' int2str(blv.npyr.final) ' nInt.f=' int2str(blv.nint.final)];[' nPyr.ccg=' int2str(blv.npyr.ccg) ' nInt.ccg=' int2str(blv.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VEn + DEn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data for structure
vdEn.stats=poissonripplemod(ismember(poissonripplemod(:,1:3),[VEn;DEn],'rows'),:);
vdEn.finaltype=finaltype(ismember(poissonripplemod(:,1:3),[VEn;DEn],'rows'));
vdEn.ccgtype=ccgtype(ismember(poissonripplemod(:,1:3),[VEn;DEn],'rows'));

% counts
vdEn.nall.final=size(vdEn.stats,1);
vdEn.nall.ccg=sum(vdEn.ccgtype~=0);
vdEn.nint.final=sum(vdEn.finaltype==2);
vdEn.nint.ccg=sum(vdEn.ccgtype==2);
vdEn.npyr.final=sum(vdEn.finaltype==1);
vdEn.npyr.ccg=sum(vdEn.ccgtype==1);

vdEn.allsig=sum(vdEn.stats(:,6)<pval/2|vdEn.stats(:,7)<pval/2);
vdEn.allsigup=sum(vdEn.stats(:,6)<pval/2);
vdEn.allsigdown=sum(vdEn.stats(:,7)<pval/2);

% int final
vdEn.intsig.final=sum((vdEn.stats(:,6)<pval/2|vdEn.stats(:,7)<pval/2)&vdEn.finaltype==2);
vdEn.intsigup.final=sum(vdEn.stats(:,6)<pval/2&vdEn.finaltype==2);
vdEn.intsigdown.final=sum(vdEn.stats(:,7)<pval/2&vdEn.finaltype==2);
%int ccg
vdEn.intsig.ccg=sum((vdEn.stats(:,6)<pval/2|vdEn.stats(:,7)<pval/2)&vdEn.ccgtype==2);
vdEn.intsigup.ccg=sum(vdEn.stats(:,6)<pval/2&vdEn.ccgtype==2);
vdEn.intsigdown.ccg=sum(vdEn.stats(:,7)<pval/2&vdEn.ccgtype==2);

% pyr final
vdEn.pyrsig.final=sum((vdEn.stats(:,6)<pval/2|vdEn.stats(:,7)<pval/2)&vdEn.finaltype==1);
vdEn.pyrsigup.final=sum(vdEn.stats(:,6)<pval/2&vdEn.finaltype==1);
vdEn.pyrsigdown.final=sum(vdEn.stats(:,7)<pval/2&vdEn.finaltype==1);
% pyr ccg
vdEn.pyrsig.ccg=sum((vdEn.stats(:,6)<pval/2|vdEn.stats(:,7)<pval/2)&vdEn.ccgtype==1);
vdEn.pyrsigup.ccg=sum(vdEn.stats(:,6)<pval/2&vdEn.ccgtype==1);
vdEn.pyrsigdown.ccg=sum(vdEn.stats(:,7)<pval/2&vdEn.ccgtype==1);

% percents
vdEn.percentall=vdEn.allsig/vdEn.nall.final*100;
vdEn.percentallup=vdEn.allsigup/vdEn.nall.final*100;
vdEn.percentalldown=vdEn.allsigdown/vdEn.nall.final*100;

vdEn.percentint.final=vdEn.intsig.final/vdEn.nint.final*100;
vdEn.percentintup.final=vdEn.intsigup.final/vdEn.nint.final*100;
vdEn.percentintdown.final=vdEn.intsigdown.final/vdEn.nint.final*100;

vdEn.percentint.ccg=vdEn.intsig.ccg/vdEn.nint.ccg*100;
vdEn.percentintup.ccg=vdEn.intsigup.ccg/vdEn.nint.ccg*100;
vdEn.percentintdown.ccg=vdEn.intsigdown.ccg/vdEn.nint.ccg*100;

vdEn.percentpyr.final=vdEn.pyrsig.final/vdEn.npyr.final*100;
vdEn.percentpyrup.final=vdEn.pyrsigup.final/vdEn.npyr.final*100;
vdEn.percentpyrdown.final=vdEn.pyrsigdown.final/vdEn.npyr.final*100;

vdEn.percentpyr.ccg=vdEn.pyrsig.ccg/vdEn.npyr.ccg*100;
vdEn.percentpyrup.ccg=vdEn.pyrsigup.ccg/vdEn.npyr.ccg*100;
vdEn.percentpyrdown.ccg=vdEn.pyrsigdown.ccg/vdEn.npyr.ccg*100;

%  figure;
subplot(2,6,11);
bar([vdEn.percentalldown*(-1) 0 vdEn.percentpyrdown.final*(-1) vdEn.percentintdown.final*(-1) 0 vdEn.percentpyrdown.ccg*(-1) vdEn.percentintdown.ccg*(-1)],'FaceColor',rgb('MidnightBlue'),'LineStyle','none');
hold on
bar([vdEn.percentallup 0 vdEn.percentpyrup.final vdEn.percentintup.final 0 vdEn.percentpyrup.ccg vdEn.percentintup.ccg],'FaceColor',rgb('Crimson'),'LineStyle','none');
set(gca,'XTickLabel',{'All' '/' 'Pyr.f' 'Int.f' '/' 'Pyr.ccg' 'Int.ccg'});
xlabel({'VEn + DEn, Poisson Ripple Mod';['nAll=' int2str(vdEn.nall.final) ' nPyr.f=' int2str(vdEn.npyr.final) ' nInt.f=' int2str(vdEn.nint.final)];[' nPyr.ccg=' int2str(vdEn.npyr.ccg) ' nInt.ccg=' int2str(vdEn.nint.ccg) ' p<' num2str(pval)]});
ylim([-40 70]);
