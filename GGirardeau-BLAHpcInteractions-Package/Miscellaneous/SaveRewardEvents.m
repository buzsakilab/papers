%SaveRewardEvents - Creates and saves reward events on the two ends of the linear track.
%
%  USAGE
%
%    SaveRewardEvents(filenameL, filenameR, rewardsLeftLFP, rewardsRightLFP)
%
%  INPUT
%    	filenameL       	file to save the rewards events on the left
%	filenameR       	file to save the rewards events on the right
%    	rewardsLeftLFP		lfp channel number for LEFT rewards
%    	rewardsRightLFP		lfp channel number for RIGHT rewards
%
%  OUTPUT
%	event files (reload w/ SetCurrentSession)
%  SEE
%
%    See also SaveEvents.
%
%	!!! Remember left and right are conventions here : it depends how you look at and plot your positions. Be careful to decide what's left and what's right in which configuration.
%
% Copyright (C) 2013 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function SaveRewardEvents (filenameL,filenameR, rewardsLeftLFP, rewardsRightLFP)


if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveRewardEvents">SaveRewardEvents</a>'' for details).');
end

rightlfp=GetLFP(rewardsRightLFP);
leftlfp=GetLFP(rewardsLeftLFP);


[rightperiods,rightin]=Threshold(rightlfp,'>',800);
[leftperiods,leftin]=Threshold(leftlfp,'>',800);

figure;
PlotXY(rightlfp,'g');
hold on
PlotXY(leftlfp,'b');
PlotIntervals([rightperiods;leftperiods],'rectangles');


right=rightperiods(:,1);
left=leftperiods(:,1);


rightevents=NewEvents(right,['Right Reward']);
leftevents=NewEvents(left,['Left Reward']);

SaveEvents(filenameL,leftevents);
SaveEvents(filenameR,rightevents);