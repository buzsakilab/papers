function BehavioralMeasuresAll

%  BehavioralMeasuresAll - Pools behavioral measures across sessions
%
%  USAGE
%
%   BehavioralMeasuresAll
%
%  OUTPUT
%
%   BehavioralMeasure.mat : variables storing data for all sessions and animals
%
%  NOTE
%
%  SEE
%  
%    BehavioralMeasures_Stats, BehavioralMeasures
%
% Gabrielle Girardeau, March 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


load('/media/Data-01/All-Rats/sessionindexing.mat');

PrePmeanVAll={};
YPrePmeanVAll={};
ratioAll={};
sessionsxml={};
ratsession=[];
rats=[];

for i=1:length(ratsessionindex)
  cd(xmlpath{i})
  xml=xmlpath{i}(end-14:end-1)
  if exist('BehavioralMeasures.mat')==2
    load('BehavioralMeasures.mat');
    PrePmeanVAll=[PrePmeanVAll;PrePmeanV];
    YPrePmeanVAll=[YPrePmeanVAll;YPrePmeanV];
    ratioAll=[ratioAll;ratio];
    sessionsxml=[sessionsxml;xml];
    rats=[rats;ratsessionindex(i,1)];
    ratsession=[ratsession;ratsessionindex(i,:)];
  end
end

cd('/media/Data-01/All-Rats/BehavioralMeasures');
save('BehavioralMeasures.mat','YPrePmeanVAll','PrePmeanVAll','ratioAll','ratsession','rats','sessionsxml','NormalizedAPpos','directions');
