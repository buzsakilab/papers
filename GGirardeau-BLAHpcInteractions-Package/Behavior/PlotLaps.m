function PlotLaps

laps=dir('*-Laps.mat')
load(laps.name);
load('airpuff.mat');
%  load('runs.mat');
pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);

%  runint=RunIntervals(runs)
load('runintervals.mat')
runint=runintervals;

runpos=Restrict(pos,runint);
figure;
PlotXY(runpos(:,1:2));
set(gcf,'Position',[355 106 1370 871]);
PlotHVLines(airpuff.loc*0.43,'h','r');
hold on


RtoLpos=Restrict(runpos,RtoLlaps);
LtoRpos=Restrict(runpos,LtoRlaps);
if ~isempty(Uturnlaps)
  Uturnpos=Restrict(runpos,Uturnlaps);
  PlotXY(Uturnpos(:,1:2),'.k');
end
PlotXY(RtoLpos(:,1:2),'.g');
PlotXY(LtoRpos(:,1:2),'.r');
