function [ratsess,ratsessxml,performance] = Performance_All

%Performance_All - Pools perfromance data across animals and plots habotuation curves across days
%
%  USAGE
%
%        [ratsess,ratsessxml,performance] = Performance_All
%
%  OUTPUT
%
%       Plots
%
%  SEE
%
%       Performance.m
%
% Gabrielle Girardeau, 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/TrainingSessionsForPerf.mat');

performance=[];

for i=1:size(xmlpath,1)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
    if exist ('Perf.mat')==2;
        load('Perf.mat');
        performance=[performance;perf];
    else
        performance=[performance;NaN]
    end
end

figure;hold on
plot(trainingday(Ratnum==8),performance(Ratnum==8),'.-','MarkerSize',14,'MarkerFaceColor',rgb('ForestGreen'),'Color',rgb('ForestGreen'))
plot(trainingday(Ratnum==9),performance(Ratnum==9),'.-','MarkerSize',14,'MarkerFaceColor',rgb('DarkGreen'),'Color',rgb('DarkGreen'))
plot(trainingday(Ratnum==10),performance(Ratnum==10),'.-','MarkerSize',14,'MarkerFaceColor',rgb('LimeGreen'),'Color',rgb('LimeGreen'))
plot(trainingday(Ratnum==11),performance(Ratnum==11),'.-','MarkerSize',14,'MarkerFaceColor',rgb('OliveDrab'),'Color',rgb('OliveDrab'))
