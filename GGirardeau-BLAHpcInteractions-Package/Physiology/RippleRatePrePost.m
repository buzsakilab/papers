function [Rate,Ratio,Gain] = RippleRatePrePost (session,presleep,postsleep)

%RippleRatePrePost - Calculates ripple rate in pre vs post sleep + Ratio btw the 2.
%
%  USAGE
%
%    [Rate,Ratio,Gain] = RippleRatePrePost (session,presleep,postsleep)
%
%    session      path to session
%    presleep     presleep epoch number
%    postsleep    postsleep epoch number
%
%  OUTPUT
%
%    Rate   .pre/.post : ripple occurrence rates
%    Ratio  pre/post ratio for ripple rate
%    Gain   Gain between pre and post rates.
%
%  October 2013, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

cd(session)
xml=session(end-13:end)
SetCurrentSession(xml,'spikes','off')
load('States.mat');

ripples=GetEvents({'Ripple peak.*'});

allPRE=RunIntervals(presleep);
allPOST=RunIntervals(postsleep);

if size(allPRE,1)>1
  allPRE=[allPRE(1,1) allPRE(end,2)];
end
if size(allPOST,1)>1
  allPOST=[allPOST(1,1) allPOST(end,2)];
end

swsPRE=Restrict(sws,allPRE);
swsPOST=Restrict(sws,allPOST);

ripplesPRE = Restrict(ripples, swsPRE,'shift','on');
ripplesPOST = Restrict (ripples, swsPOST,'shift','on');

swstimePRE=sum(swsPRE(:,2)-swsPRE(:,1))
swstimePOST=sum(swsPOST(:,2)-swsPOST(:,1))

Rate.PRE=length(ripplesPRE)/swstimePRE;
Rate.POST=length(ripplesPOST)/swstimePOST;

Ratio=(Rate.POST-Rate.PRE)/(Rate.POST+Rate.PRE);
Gain=(Rate.POST/Rate.PRE);

save([xml '-RipplesPrePost.mat'],'Rate','Ratio','Gain');

cd('/media/Data-01/All-Rats/AllRats-RippleRates');
figure('Position',[817 516 865 371]);
subplot(1,4,3)
bar([Rate.PRE Rate.POST],'LineStyle','none');
ylabel('Mean Ripple Rate');
subplot(1,4,1:2)
plot(ripplesPRE,1:length(ripplesPRE),'b');
hold on;
plot(ripplesPOST,1:length(ripplesPOST),'r');
subplot(1,4,4)
bar([swstimePRE swstimePOST],'LineStyle','none');
ylabel('SWS time')
suptitle([xml ' Ripple ratio : ' num2str(Ratio)]);
plot2svg([xml '-SWSRipples.svg'],gcf);

%  figure;
%  plot(ripples,cumsum(ones(length(ripples),1)));
%  PlotColorIntervals(swsPRE,'Green','v');
%  PlotColorIntervals(swsPOST,'Yellow','v');

