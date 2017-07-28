function RippleRatePrePost_All

%RippleRatePrePost_All - Aggregates ripple rate data across rats and sessions
%
%  USAGE
%
%    RippleRatePrePost_All
%
%  OUTPUT
%
%    Saved variable
%
%  October 2013, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/sessionindexing.mat');

rate.pre=[];
rate.post=[];

ratsess=[];

for i=1:size(ratsessionindex)
    currentsession=xmlpath{i}
    cd(currentsession);
    xml=currentsession(end-14:end-1);
    if exist ([xml '-RipplesPrePost.mat'])==2;
        load([xml '-RipplesPrePost.mat']);
        rate.pre=[rate.pre;Rate.PRE];
        rate.post=[rate.post;Rate.POST];
        ratsess=[ratsess;ratsessionindex(i,:)];
    end
end

figure;
boxplot([rate.pre rate.post]);
hold on;
for i=1:length(rate.pre)
    plot([1 2],[rate.pre(i) rate.post(i)],'k');
end

[h,p.stats]=signrank(rate.pre,rate.post,'tail','right')