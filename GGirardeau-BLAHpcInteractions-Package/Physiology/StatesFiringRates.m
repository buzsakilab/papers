function StatesFR = StatesFiringRates (session)

%  StatesFiringRates - Calculates and plots firing rates for different States (Wake, REM and SWS) for one session
%
%  USAGE
%
%    	StatesFR = StatesFiringRates(session)
%
%  INPUT
%
%	session	      path to session
%
%  OUTPUT
%
%  	StatesFR 	Matrix containing : RAt / Session / ShankNumber / CellNumber / FR-REM / FR-SWS / FR-Wake / FR-Drowsy / FR-All
%
%
%
% Copyright (C) 2013-14 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

StatesFR=[];

cd(session);
xml=session(end-13:end);
SetCurrentSession(xml);

load('/media/Data-01/All-Rats/sessionindexing.mat');
ratsess = ratsessionindex(strcmp(xmlpath,[session '/']),:)

% Get All Spikes
spikes1=GetSpikeTimes('output','numbered');
spikes2=GetSpikeTimes('output','full');
spikes=[spikes2 spikes1(:,2)];
spikes(spikes(:,3)==0,:)=[];
spikes(spikes(:,3)==1,:)=[];

% Remnumber IDs after removing MUA and artifacts
ids=unique(spikes(:,4));
for ii=1:length(ids)
  spikes(spikes(:,4)==ids(ii),4)=ii;
end
idx=unique(spikes(:,2:4),'rows'); % shank / cell/ ID

if min(idx(:,3))==1 & sum(diff(idx(:,3))-1)==0
  disp('Renumbering OK')
else
  error('Problem with index : aborting')
end

% Brain states
load('States.mat');
remintervals=Rem;
swsintervals=sws;
explointervals=wake;
drowsyintervals=drowsy;

for id=1:size(idx,1)
    unittimestamps=spikes(spikes(:,4)==id);
    % rem sleep
    [statusrem,interval,~]=InIntervals(unittimestamps,remintervals);
    remspikes=unittimestamps(statusrem);
    remfreq=length(remspikes)/sum(remintervals(:,2)-remintervals(:,1));
    % swsintervals
    [statussws,interval,~]=InIntervals(unittimestamps,swsintervals);
    swsspikes=unittimestamps(statussws);
    swsfreq=length(swsspikes)/sum(swsintervals(:,2)-swsintervals(:,1));
    % wake
    [statusexplo,interval,~]=InIntervals(unittimestamps,explointervals);
    explospikes=unittimestamps(statusexplo);
    explofreq=length(explospikes)/sum(explointervals(:,2)-explointervals(:,1));
    % drowsy
    [statusdrowsy,interval,~]=InIntervals(unittimestamps,drowsyintervals);
    drowsyspikes=unittimestamps(statusdrowsy);
    drowsyfreq=length(drowsyspikes)/sum(drowsyintervals(:,2)-drowsyintervals(:,1));
    % all states
    start=min([remintervals(:,1);swsintervals(:,1);explointervals(:,1);drowsyintervals(:,1)]);
    stop=max([remintervals(:,2);swsintervals(:,2);explointervals(:,2);drowsyintervals(:,2)]);
    allstatesfreq=length(unittimestamps)/(stop-start);
    % fill in the output
    StatesFR(id,:)=[ratsess idx(id,:) remfreq swsfreq explofreq drowsyfreq allstatesfreq];
end

