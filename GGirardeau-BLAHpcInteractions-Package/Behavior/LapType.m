function [safelaps,aplaps]=LapType(session)
% [safelaps,airpufflaps]=LapType(LtoRLaps,RtoLlaps,airpuffdir,runintervals)
% classifies safelaps and airpuff laps using airpuffdirection - Should be 3 runintervals (preun-run-postrun)
% OUTPUT : safelaps/airpufflaps		.all	.prerun	.run .postrun

cd(session);
xml=session(end-13:end);

load([xml '-Laps.mat']);
load('States.mat');
load('runintervals.mat');
load('airpuff.mat')
git st
if size(runintervals,1)<3
  warning('Should be 3 runintervals : wrong size')
end

if strcmp(airpuff.dir,'LtoR')
  aplaps.all = LtoRlaps ;
  aplaps.prerun = Restrict(LtoRlaps,runintervals(1,:));
  aplaps.run = Restrict (LtoRlaps,runintervals(2,:));
  aplaps.postrun = Restrict (LtoRlaps,runintervals(3,:));

  safelaps.all = RtoLlaps ;
  safelaps.prerun = Restrict(RtoLlaps,runintervals(1,:));
  safelaps.run = Restrict (RtoLlaps,runintervals(2,:));
  safelaps.postrun = Restrict (RtoLlaps,runintervals(3,:));

elseif strcmp(airpuff.dir,'RtoL')
  aplaps.all = RtoLlaps ;
  aplaps.prerun = Restrict(RtoLlaps,runintervals(1,:));
  aplaps.run = Restrict (RtoLlaps,runintervals(2,:));
  aplaps.postrun = Restrict (RtoLlaps,runintervals(3,:));

  safelaps.all = LtoRlaps ;
  safelaps.prerun = Restrict(LtoRlaps,runintervals(1,:));
  safelaps.run = Restrict (LtoRlaps,runintervals(2,:));
  safelaps.postrun = Restrict (LtoRlaps,runintervals(3,:));
end

save([xml '-LapType.mat'],'safelaps','aplaps');