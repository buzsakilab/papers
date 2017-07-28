function [NormalizedAPpos,xmlIndex,directions] = AirpuffPositionDistribution

%  AirpuffPositionDistribution - Normalizes airpuff position across sessions 
%
%  USAGE
%
%    [NormalizedAPpos,xmlIndex,directions] = AirpuffPositionDistribution
%
%  OUTPUT
%
%    NormalizedAPpos    Normalized airpuff position   
%    xmlIndex           Corresponding list of sessions  
%    directions         Direction in which the airpuff is given
%
%  NOTE
%  
%  SEE
%
%    See also : AirpuffPositionDistribution_Plot, ZeroToOne (FMAToolbox)
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


load('/media/Data-01/All-Rats/sessionindexing.mat');


NormalizedAPpos=[];
xmlIndex=[];
ratsess=[];
directions=[]; %LtoR = 0, RtoL=1
for i=1:size(ratsessionindex)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  if exist ([xml '-Laps.mat'])==2 & exist('airpuff.mat')==2
    load([xml '-Laps.mat']);
    load('airpuff.mat');
    ap=airpuff.loc*0.43;
    
    [~,normap]=ZeroToOne([leftLimit rightLimit],ap);
    NormalizedAPpos=[NormalizedAPpos;normap];
    xmlIndex=[xmlIndex;xml];
    if strcmp(airpuff.dir,'LtoR');
        directions=[directions;0];
    elseif strcmp(airpuff.dir,'RtoL');
        directions=[directions;1];
    end
    ratsess=[ratsess;ratsessionindex(i,:)];
  end
end

cd('/media/Data-01/All-Rats/BehavioralMeasures');
save('NormAPpos.mat','NormalizedAPpos','xmlIndex','directions','ratsess');
