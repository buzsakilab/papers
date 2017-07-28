function [laps,runtime,perf,xml] = Performance(session,run)

%Performance - Calculates the Performance for a single session (number of runned laps/time)
%
%  USAGE
%
%    [laps,runtime,perf,xml] = Performance(session,run)
%
%    run            run subsession number
%    session        session path
%
%  OUTPUT
%
%    laps           number of laps
%    runtime        total runtime
%    perf           number of laps/runtim*100
%
%  SEE
%
%    Performance_All
%
% Gabrielle Girardeau, 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

cd(session);
xml=session(end-13:end);
SetCurrentSession('filename',xml,'spikes','off');

runInterval = RunIntervals(run);
load([xml '-Laps.mat']);
runlrw=Restrict(LtoRlaps(:,1),runInterval);

laps=size(runlrw,1);
runtime=runInterval(2)-runInterval(1);
perf=(laps/runtime)*100;

save('Perf.mat','perf','laps','runtime');