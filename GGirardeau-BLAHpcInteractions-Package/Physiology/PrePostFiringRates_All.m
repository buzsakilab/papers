function PrePostFRAll = PrePostFiringRates_All
%PrePostFiringRatesAll - Aggregates firing rates of individual cells over rats and sessions for each of the three states, before and after learning.
%
%  USAGE
%  
%    function PrePostFRAll = PrePostFiringRates_All
%
%  OUTPUT
%
%    PrePostFRAll : matrix [Rat/session/shank/unit/id/preremrate/postremrate/preswsrate/postswsrate/prewakerate/postwakerate/prerunrate/runrate/postrunrate]
%
%  NOTES
%
%  SEE
%  
%    PrePostFR
%
% 2016 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.



load('/media/Data-01/All-Rats/sessionindexing.mat');


PrePostFRAll=[];

for i=1:size(ratsessionindex)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  if exist ([xml '-PrePostFR.mat'])==2;
    load([xml '-PrePostFR.mat']);
    PrePostFRAll=[PrePostFRAll;PrePostFR];
  end
end

cd('/media/Data-01/All-Rats')
save('AllRats-PrePostFRAll.mat','PrePostFRAll');