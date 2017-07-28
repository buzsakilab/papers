function VelocityMaps_All(varargin)

%VelocityMaps_All - Calculates average speed curves in the safe and airpuff directions across rats/sessions
%
%  USAGE
%
%    [safemap,apmap,airP] = VelocityMaps(session)
%
%    session    path to session
%
%  OUTPUT
%  
%   safemap         speed curves for safe runs (.prerun .run .postrun)
%   apmap           speed curves for airpuff runs
%   airP            normalized airpuff position
%
%  SEE
%
%    VelocityMaps_All
%
% Gabrielle Girardeau, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


%  Defaults
postlaps = 'all' ; % or 'firsts' for first 5 laps.

% Check varargin
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters  (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'postlaps',
			postlaps = varargin{i+1};
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end


load('/media/Data-01/All-Rats/sessionindexing.mat');

alignedMatrix.run.ap=[];
alignedMatrix.run.safe=[];
alignedMatrix.prerun.ap=[];
alignedMatrix.prerun.safe=[];
alignedMatrix.postrun.ap=[];
alignedMatrix.postrun.safe=[];

meanvel.run.safe=[];
meanvel.run.ap=[];
meanvel.prerun.safe=[];
meanvel.prerun.ap=[];
meanvel.postrun.safe=[];
meanvel.postrun.ap=[];

allxml=[];
allratsess=[];

for i=1:size(ratsessionindex)
    currentsession=xmlpath{i}
    cd(currentsession);
    xml=currentsession(end-14:end-1);
    if exist([xml '-VelocityCurves.mat'])==2
        load([xml '-VelocityCurves.mat']);
        load('airpuff.mat');

        % matrix of normalized x positions from -1 to 1 and corresponding NaN matrix to fill with values
        xbasematrix = [flip(apmap.run.x(2:50)')*(-1)';apmap.run.x'];
        xbasematrix = round(xbasematrix*10000)./10000;
        NaNbase=NaN(1,99);
        
        
        x=apmap.run.x;
        [~,closestx]=FindClosest(x,airP);
        APcenteredx=x-closestx;
        APcenteredx=round(APcenteredx'*10000)./10000;
        
        velcurve.run.ap=NaNbase;
        velcurve.run.safe=NaNbase;
        velcurve.prerun.ap=NaNbase;
        velcurve.prerun.safe=NaNbase;
        velcurve.postrun.ap=NaNbase;
        velcurve.postrun.safe=NaNbase;
        
        meanvel.run.ap=[meanvel.run.ap;mean(apmap.run.z)];
        meanvel.run.safe=[meanvel.run.safe;mean(safemap.run.z)];
        meanvel.prerun.ap=[meanvel.prerun.ap;mean(apmap.prerun.z)];
        meanvel.prerun.safe=[meanvel.prerun.safe;mean(safemap.prerun.z)];
        if strcmp(postlaps,'all')
            meanvel.postrun.ap=[meanvel.postrun.ap;mean(apmap.postrun.all.z)];
            meanvel.postrun.safe=[meanvel.postrun.safe;mean(safemap.postrun.all.z)];
        elseif strcmp(postlaps,'firsts')
            meanvel.postrun.ap=[meanvel.postrun.ap;mean(apmap.postrun.firsts.z)];
            meanvel.postrun.safe=[meanvel.postrun.safe;mean(safemap.postrun.firsts.z)];
        end    

        if strcmp(airpuff.dir,'RtoL')	%flip airpuff traj
            flipx.run.ap=flip(APcenteredx*(-1));
            velcurve.prerun.ap(ismember(xbasematrix,flipx.run.ap))=flip(apmap.prerun.z);
            velcurve.prerun.safe(ismember(xbasematrix,APcenteredx))=safemap.prerun.z;
            velcurve.run.ap(ismember(xbasematrix,flipx.run.ap))=flip(apmap.run.z);
            velcurve.run.safe(ismember(xbasematrix,APcenteredx))=safemap.run.z;
            if strcmp(postlaps,'all')
                if ~isempty(apmap.postrun.all.z)
                    velcurve.postrun.ap(ismember(xbasematrix,flipx.run.ap))=flip(apmap.postrun.all.z);
                end
                velcurve.postrun.safe(ismember(xbasematrix,APcenteredx))=safemap.postrun.all.z;
            elseif strcmp(postlaps,'firsts')    
                if ~isempty(apmap.postrun.firsts.z)
                    velcurve.postrun.ap(ismember(xbasematrix,flipx.run.ap))=flip(apmap.postrun.firsts.z);
                end
                velcurve.postrun.safe(ismember(xbasematrix,APcenteredx))=safemap.postrun.firsts.z;
            end
        else
            flipx.run.safe=flip(APcenteredx*(-1));
            if ~isempty(safemap.prerun.x)
                velcurve.prerun.safe(ismember(xbasematrix,APcenteredx))=flip(safemap.prerun.z);
                velcurve.prerun.ap(ismember(xbasematrix,APcenteredx))=apmap.prerun.z;
            end
            velcurve.run.safe(ismember(xbasematrix,APcenteredx))=flip(safemap.run.z);
            velcurve.run.ap(ismember(xbasematrix,APcenteredx))=apmap.run.z;
            if strcmp(postlaps,'all')
                velcurve.postrun.safe(ismember(xbasematrix,APcenteredx))=flip(safemap.postrun.all.z);
                velcurve.postrun.ap(ismember(xbasematrix,APcenteredx))=apmap.postrun.all.z;
            elseif strcmp(postlaps,'firsts')    
                velcurve.postrun.safe(ismember(xbasematrix,APcenteredx))=flip(safemap.postrun.firsts.z);
                velcurve.postrun.ap(ismember(xbasematrix,APcenteredx))=apmap.postrun.firsts.z;
            end
        end
        alignedMatrix.run.ap=[alignedMatrix.run.ap;velcurve.run.ap];
        alignedMatrix.run.safe=[alignedMatrix.run.safe;velcurve.run.safe];
        alignedMatrix.prerun.ap=[alignedMatrix.prerun.ap;velcurve.prerun.ap];
        alignedMatrix.prerun.safe=[alignedMatrix.prerun.safe;velcurve.prerun.safe];
        alignedMatrix.postrun.ap=[alignedMatrix.postrun.ap;velcurve.postrun.ap];
        alignedMatrix.postrun.safe=[alignedMatrix.postrun.safe;velcurve.postrun.safe];
        
        allxml=[allxml;xml];
    end
end

figure('Position',[451 623 1348 270]);
subplot(3,1,2);hold on;
plot(xbasematrix,nanmean(alignedMatrix.run.ap),'r');
plot(xbasematrix,nanmean(alignedMatrix.run.ap)+nansem(alignedMatrix.run.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.run.ap)-nansem(alignedMatrix.run.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.run.safe));
plot(xbasematrix,nanmean(alignedMatrix.run.safe)+nansem(alignedMatrix.run.safe),'b:');
plot(xbasematrix,nanmean(alignedMatrix.run.safe)-nansem(alignedMatrix.run.safe),'b:');
ylim([10 90]);
PlotHVLines(0,'v','r');
%  subplot(3,2,4)
%  plot(xbasematrix,nanmean(alignedMatrix.run.safe));
%  ylim([10 90]);
%  PlotHVLines(0,'v','b');

subplot(3,1,1);hold on;
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap),'r');
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap)+nansem(alignedMatrix.prerun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap)-nansem(alignedMatrix.prerun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.prerun.safe));
plot(xbasematrix,nanmean(alignedMatrix.prerun.safe)+nansem(alignedMatrix.prerun.safe),'b:');
plot(xbasematrix,nanmean(alignedMatrix.prerun.safe)-nansem(alignedMatrix.prerun.safe),'b:');
ylim([10 90]);
PlotHVLines(0,'v','r');
%  subplot(3,2,2)
%  plot(xbasematrix,nanmean(alignedMatrix.prerun.safe));
%  ylim([10 90]);
%  PlotHVLines(0,'v','b');

subplot(3,1,3);hold on;
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap),'r');
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap)+nansem(alignedMatrix.postrun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap)-nansem(alignedMatrix.postrun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.safe)+nansem(alignedMatrix.postrun.safe),'b:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.safe)-nansem(alignedMatrix.postrun.safe),'b:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.safe));
ylim([10 90]);
PlotHVLines(0,'v','r');
%  subplot(3,2,6)
%  plot(xbasematrix,nanmean(alignedMatrix.postrun.safe));
%  ylim([10 90]);
%  PlotHVLines(0,'v','b');

figure;
subplot(1,3,1)
bar([mean(meanvel.prerun.ap) mean(meanvel.prerun.safe)]);
subplot(1,3,2)
bar([mean(meanvel.run.ap) mean(meanvel.run.safe)]);
subplot(1,3,3)
bar([nanmean(meanvel.postrun.ap) nanmean(meanvel.postrun.safe)]);

DZmeanvel.prerun.ap=mean(alignedMatrix.prerun.ap(:,44:50),2);
DZmeanvel.prerun.safe=mean(alignedMatrix.prerun.safe(:,44:50),2);
DZmeanvel.run.ap=mean(alignedMatrix.run.ap(:,44:50),2);
DZmeanvel.run.safe=mean(alignedMatrix.run.safe(:,44:50),2);
DZmeanvel.postrun.ap=mean(alignedMatrix.postrun.ap(:,44:50),2);
DZmeanvel.postrun.safe=mean(alignedMatrix.postrun.safe(:,44:50),2);

figure;
bar([nanmean(DZmeanvel.prerun.safe) nanmean(DZmeanvel.prerun.ap) nanmean(DZmeanvel.run.safe) nanmean(DZmeanvel.run.ap) nanmean(DZmeanvel.postrun.safe) nanmean(DZmeanvel.postrun.ap)])

[h,p,stats]=ranksum(DZmeanvel.prerun.ap,DZmeanvel.postrun.ap);

[h,p]=ranksum(DZmeanvel.prerun.safe-DZmeanvel.prerun.ap,DZmeanvel.postrun.safe-DZmeanvel.postrun.ap);
[h,p]=ranksum(DZmeanvel.prerun.safe,DZmeanvel.prerun.ap)
[h,p]=ranksum(DZmeanvel.postrun.safe,DZmeanvel.postrun.ap)

figure;hold on;
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap),'b');
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap)+nansem(alignedMatrix.prerun.ap),'b:');
plot(xbasematrix,nanmean(alignedMatrix.prerun.ap)-nansem(alignedMatrix.prerun.ap),'b:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap),'r');
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap)+nansem(alignedMatrix.postrun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.postrun.ap)-nansem(alignedMatrix.postrun.ap),'r:');
plot(xbasematrix,nanmean(alignedMatrix.run.ap),'k');
plot(xbasematrix,nanmean(alignedMatrix.run.ap)+nansem(alignedMatrix.run.ap),'k:');
plot(xbasematrix,nanmean(alignedMatrix.run.ap)-nansem(alignedMatrix.run.ap),'k:');

keyboard

for i=1:53
    if i==1
        M.prerun.ap(i,:)=NaN(1,50);
        M.prerun.safe(i,:)=NaN(1,50);
    else
        M.prerun.safe(i,:)=alignedMatrix.prerun.safe(i,~isnan(alignedMatrix.prerun.safe(i,:)));
        M.prerun.ap(i,:)=alignedMatrix.prerun.ap(i,~isnan(alignedMatrix.prerun.ap(i,:)));
    end
    M.run.safe(i,:)=alignedMatrix.run.safe(i,~isnan(alignedMatrix.run.safe(i,:)));
    M.run.ap(i,:)=alignedMatrix.run.ap(i,~isnan(alignedMatrix.run.ap(i,:)));
    M.postrun.safe(i,:)=alignedMatrix.postrun.safe(i,~isnan(alignedMatrix.postrun.safe(i,:)));
    if i==3
        M.postrun.ap(i,:)=NaN(1,50);
    else
        M.postrun.ap(i,:)=alignedMatrix.postrun.ap(i,~isnan(alignedMatrix.postrun.ap(i,:)));
    end
end

allcurves.pre=[M.prerun.safe;M.prerun.ap];
allcurves.post=[M.postrun.safe;M.postrun.ap];
allcurves.run=[M.run.safe;M.run.ap];

figure;
hold on;
subplot(1,3,1)
OrderedImagesc(M.prerun.ap,'OrderType','none','newfig','off')
colorbar
subplot(1,3,2)
OrderedImagesc(M.run.ap,'OrderType','none','newfig','off')
colorbar
subplot(1,3,3)
OrderedImagesc(M.postrun.ap,'OrderType','none','newfig','off')
colorbar

figure;
hold on;
subplot(1,3,1)
OrderedImagesc(alignedMatrix.prerun.ap,'OrderType','none','newfig','off')
colorbar
subplot(1,3,2)
OrderedImagesc(alignedMatrix.run.ap,'OrderType','none','newfig','off')
colorbar
subplot(1,3,3)
OrderedImagesc(alignedMatrix.postrun.ap,'OrderType','none','newfig','off')
colorbar

%  figure;hold on;
%  subplot(1,3,1);hold on;
%  plot(nanmean(allcurves.pre),'b');
%  plot(nanmean(allcurves.pre)+nansem(allcurves.pre),'b:');
%  plot(nanmean(allcurves.pre)-nansem(allcurves.pre),'b:');
%  plot(nanmean(allcurves.post),'r');
%  plot(nanmean(allcurves.post)+nansem(allcurves.post),'r:');
%  plot(nanmean(allcurves.post)-nansem(allcurves.post),'r:');
%  plot(nanmean(allcurves.run),'k');
%  plot(nanmean(allcurves.run)+nansem(allcurves.run),'k:');
%  plot(nanmean(allcurves.run)-nansem(allcurves.run),'k:');
%  ylim([0 80])
%  xlim([0 51])
%  xlabel('AllCurves - PRe/run/post')
%  subplot(1,3,2);hold on;
%  plot(nanmean(M.prerun.safe),'b');
%  plot(nanmean(M.prerun.safe)+nansem(M.prerun.safe),'b:');
%  plot(nanmean(M.prerun.safe)-nansem(M.prerun.safe),'b:');
%  plot(nanmean(M.postrun.safe),'r');
%  plot(nanmean(M.postrun.safe)+nansem(M.postrun.safe),'r:');
%  plot(nanmean(M.postrun.safe)-nansem(M.postrun.safe),'r:');
%  plot(nanmean(M.run.safe),'k');
%  plot(nanmean(M.run.safe)+nansem(M.run.safe),'k:');
%  plot(nanmean(M.run.safe)-nansem(M.run.safe),'k:');
%  ylim([0 80])
%  xlim([0 51])
%  xlabel('Safe Runs - Pre/Run/POst')
%  subplot(1,3,3);hold on;
%  plot(nanmean(M.prerun.ap),'b');
%  plot(nanmean(M.prerun.ap)+nansem(M.prerun.ap),'b:');
%  plot(nanmean(M.prerun.ap)-nansem(M.prerun.ap),'b:');
%  plot(nanmean(M.postrun.ap),'r');
%  plot(nanmean(M.postrun.ap)+nansem(M.postrun.ap),'r:');
%  plot(nanmean(M.postrun.ap)-nansem(M.postrun.ap),'r:');
%  plot(nanmean(M.run.ap),'k');
%  plot(nanmean(M.run.ap)+nansem(M.run.ap),'k:');
%  plot(nanmean(M.run.ap)-nansem(M.run.ap),'k:');
%  ylim([0 80])
%  xlim([0 51])
%  xlabel('Airpuff runs - Pre/Run/Post')