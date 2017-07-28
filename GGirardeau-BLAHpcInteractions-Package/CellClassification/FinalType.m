function finalType = FinalType (SpikeParameters,kmeanstype);

%FinalType - Creates a variable for cell types using kmeans clustering results or CCG monosynaptic classification when available.
%
%  USAGE
%
%    finaltype = FinalType (SpikeParameters,kmeanstype);
%
%    SpikeParameters      A matrix Rat / Session / Shank / Cell / TtoP /  invF / invF2 / RatioTtoP / asym / meanFR / CCGType
%
%  OUTPUT
%
%    finaltype           A matrix Rat / sessions/ shank /unit / final type : 1=pyramidal cell, 2= interneuron.
%
%  NOTE
%
%  SEE
%
%    See also : SpikesParameters, SpikeParametersAll, CellClassifkmeans
%
% Gabrielle Girardeau, July 2014
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

finalType = kmeanstype;
%  Replace kmeans classification results by CCG monosynaptic classification wheh available.
finalType(SpikeParameters(:,end)~=0,5)=SpikeParameters(SpikeParameters(:,end)~=0,end);
save('/meida/Data-01/All-Rats/AllRats-FinalType.mat','finalType');

