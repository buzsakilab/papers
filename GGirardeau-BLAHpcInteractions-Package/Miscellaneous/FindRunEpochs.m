function FindRunEpochs (session,runs)
%Finds periods of high and low activity on the track and homecage

cd(session);
xml=session(end-13:end);
SetCurrentSession(xml,'spikes','off');
load('States.mat');

pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);
vel=LinearVelocity(pos,10);

if runs==0
  runint=[];
else
  runint=RunIntervals(runs);
end

[runtimes,in]=Threshold(vel,'>=',5,'min',2,'max',0.5);
[quiettimes,in]=Threshold(vel,'<=',3,'min',2,'max',0.5);

if ~isempty(runint)
  trackquiettimes=Restrict(quiettimes,runint);
  trackruntimes=Restrict(runtimes,runint);
  cagequiettimes=ExcludeIntervals(quiettimes,runint);
  cagequiettimes=ExcludeIntervals(cagequiettimes,sws);
  cagequiettimes=ExcludeIntervals(cagequiettimes,Rem);

  save([xml '-TrackRunTimes.mat'],'trackruntimes');
  save([xml '-QuietTimes.mat'],'trackquiettimes','cagequiettimes');
else
  cagequiettimes=ExcludeIntervals(quiettimes,sws);
  cagequiettimes=ExcludeIntervals(cagequiettimes,Rem);
  save([xml '-QuietTimes.mat'],'cagequiettimes');
end

figure;
PlotXY(pos(:,1:2));
hold on;
PlotColorIntervals(trackruntimes,'Grey','v')
