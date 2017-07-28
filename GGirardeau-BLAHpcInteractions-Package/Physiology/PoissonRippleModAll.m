function allmod = PoissonRippleModAll
% PoissonRippleModAll - Aggregates Poisson matlab variables for all rats
%
%  USAGE
%
%    allmod = PoissonRippleModAll %
%
%  OUTPUT
%
%    allmod : matrix [Rat / session / Shank / Cell / ID / pIncrease / pDecrease / Surprise]
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

load('/media/Data-01/All-Rats/sessionindexing.mat')

cd (session)
xml=session(end-13:end)
SetCurrentSession(xml,'spikes','off');

ratsess = ratsessionindex(strcmp(xmlpath,[session '/']),:)

if exist([xml '-PoissonRippleMod.mat'])
  load([xml '-PoissonRippleMod.mat'])
end
ratsess=repmat(ratsess,size(poissonripplemod,1),1);

allmod=[ratsess poissonripplemod];
ripples=GetRippleEvents;