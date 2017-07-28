function CheckLaps(session)

%  CheckLaps - Quick checks tha accurracy of lap classification (Left to Right/Right to Left/Uturn) in the stored variable Laps.m
%
%  USAGE
%
%  CheckLaps (session)
%
%  OUTPUT
%
%   Plot
%
%  NOTE
%
%  SEE
%  
%    ClassifyLaps
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

cd(session);
xml=session(end-13:end);

load([xml '-Laps.mat']);
figure;

PlotColorIntervals(LtoRlaps,'Green','v');
PlotColorIntervals(RtoLlaps,'Orange','v');
if ~isempty(Uturnlaps)
  PlotColorIntervals(Uturnlaps,'Grey','v');
end
xlabel(xml);
