%  SaveAirpuffEvents - Saves airpuff events to an event file.
%
%  INPUT
%
%     filename	name of the event file (typically Rat-date.puf.evt)
%     airpufflfp	channel number to be used to detect airpuff events


function airpuffevt = SaveAirpuffEvents (filename, airpufflfp)

lfp=GetLFP(airpufflfp);
[periods, in]=Threshold(lfp,'>',3000,'max',0.5);
airpuffevt=NewEvents(periods(:,1),'Airpuff')

figure;
PlotXY(lfp);
hold on
PlotColorIntervals(periods,'Red','v');
plot(periods(:,1),ones(size(periods,1),1)*200,'*r')
PlotHVLines(3000,'h','k')


SaveEvents(filename,airpuffevt);

