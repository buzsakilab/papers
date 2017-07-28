function kmeanstype = CellClassifkmeans (SpikeParameters)

%PlotSpikeParameters - Plots Spike Parameters
%
%  USAGE
%
%    PlotSpikeParameters (AllParamaters,structure(s))
%
%    SpikeParameters      A matrix Rat / Session / Shank / Cell / TtoP /  invF / invF2 / RatioTtoP / asym / meanFR / CCGType
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'newfig'		'on' plots in a new figure, 'off' plots in existing figure (Default 'on')
%     'rats'		'all' or 'single'
%    =========================================================================
%
%  OUTPUT
%
%    Figure(s)
%
%
%  NOTE
%
%  SEE
%
%    See also : SpikesParameters, SpikeParametersAll.
%
% Gabrielle Girardeau, July 2014
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
structure=[];
newfig='on';
clusters=[];
rats='all';

% Check number of inputs
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

%  % Check varargin
%  if mod(length(varargin),2) ~= 0,
%    error('Incorrect number of parameters  (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
%  end

%  % Parse options
%  for i = 1:2:length(varargin),
%  	if ~ischar(varargin{i}),
%  		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
%  	end
%  	switch(lower(varargin{i})),
%  		case 'structure',
%  			structure = varargin{i+1};
%  		case 'newfig'
%  			newfig = varargin{i+1};
%  		case 'clusters'
%  			clusters = varargin{i+1};
%  		case 'rats'
%  			rats = varargin{i+1};
%  		otherwise,
%  			error(['Unknown property ' num2str(varargin{i})]);
%  	end
%  end


TtoP=5;
invF=6; % Adrien's : do not use
invF2=7; % Eran's : see Resonance paper.
asym=9;
meanFR=10;
type=11;

% load structures
load('/media/Data-01/All-Rats/Structures/structures.mat');

% group structures HARD CODING
basolat=[BLA;BLV;LaDL;BMP;BMA];
central=CeCM;
olfact=[Pir;VEn;DEn];
hippo=Hpc;

% initialize variable
kmeanstype=[SpikeParameters(:,1:4) zeros(size(SpikeParameters,1),1)];

SpikeP.basolat=SpikeParameters(ismember(SpikeParameters(:,1:3),basolat,'rows'),:);
SpikeP.central=SpikeParameters(ismember(SpikeParameters(:,1:3),central,'rows'),:);
SpikeP.olfact=SpikeParameters(ismember(SpikeParameters(:,1:3),olfact,'rows'),:);
SpikeP.hippo=SpikeParameters(ismember(SpikeParameters(:,1:3),hippo,'rows'),:);

km.basolat=kmeans(SpikeP.basolat(:,[TtoP invF2]),2);
if sum(km.basolat==2)>sum(km.basolat==1)
  km.basolat=km.basolat+2;km.basolat(km.basolat==4)=1;km.basolat(km.basolat==3)=2;
end
km.central=kmeans(SpikeP.central(:,[TtoP invF2]),2);
if sum(km.central==2)>sum(km.central==1)
  km.central=km.central+2;km.central(km.central==4)=1;km.central(km.central==3)=2;
end
km.olfact=kmeans(SpikeP.olfact(:,[TtoP invF2]),2);
if sum(km.olfact==2)>sum(km.olfact==1)
  km.olfact=km.olfact+2;km.olfact(km.olfact==4)=1;km.olfact(km.olfact==3)=2;
end
km.hippo=kmeans(SpikeP.hippo(:,[TtoP invF2]),2);
if sum(km.hippo==2)>sum(km.hippo==1)
  km.hippo=km.hippo+2;km.hippo(km.hippo==4)=1;km.hippo(km.hippo==3)=2;
end

kmeanstype(ismember(SpikeParameters(:,1:3),basolat,'rows'),5)=km.basolat;
kmeanstype(ismember(SpikeParameters(:,1:3),central,'rows'),5)=km.central;
kmeanstype(ismember(SpikeParameters(:,1:3),olfact,'rows'),5)=km.olfact;
kmeanstype(ismember(SpikeParameters(:,1:3),hippo,'rows'),5)=km.hippo;


%  figure;
%  scatter3(SpikeP.basolat(km.basolat==2,TtoP),SpikeP.basolat(km.basolat==2,invF2),SpikeP.basolat(km.basolat==2,meanFR),300,rgb('Tomato'),'.');
%  hold on
%  scatter3(SpikeP.basolat(km.basolat==1,TtoP),SpikeP.basolat(km.basolat==1,invF2),SpikeP.basolat(km.basolat==1,meanFR),300,rgb('RoyalBlue'),'.');

figure;
scatter(SpikeP.basolat(km.basolat==1,TtoP),SpikeP.basolat(km.basolat==1,invF2),300,rgb('Tomato'),'.');
hold on
scatter(SpikeP.basolat(km.basolat==2,TtoP),SpikeP.basolat(km.basolat==2,invF2),300,rgb('RoyalBlue'),'.');
xlabel('basolat');

figure;
scatter(SpikeP.central(km.central==1,TtoP),SpikeP.central(km.central==1,invF2),300,rgb('Tomato'),'.');
hold on
scatter(SpikeP.central(km.central==2,TtoP),SpikeP.central(km.central==2,invF2),300,rgb('RoyalBlue'),'.');
xlabel('central');

figure;
scatter(SpikeP.olfact(km.olfact==1,TtoP),SpikeP.olfact(km.olfact==1,invF2),300,rgb('Tomato'),'.');
hold on
scatter(SpikeP.olfact(km.olfact==2,TtoP),SpikeP.olfact(km.olfact==2,invF2),300,rgb('RoyalBlue'),'.');
xlabel('olfact');

figure;
scatter(SpikeP.hippo(km.hippo==1,TtoP),SpikeP.hippo(km.hippo==1,invF2),300,rgb('Tomato'),'.');
hold on
scatter(SpikeP.hippo(km.hippo==2,TtoP),SpikeP.hippo(km.hippo==2,invF2),300,rgb('RoyalBlue'),'.');
xlabel('hippo');

%  %%% TtoP vs SpikeWidth (invF2)
%  if strcmp(newfig,'on')
%    figure;
%  end
%  subplot(3,3,[2 3 5 6]);
%    hold on;
%    scatter(SpikeParameters(SpikeParameters(:,type)==0,TtoP)+(rand(sizeUndef,1)-0.5)./50,SpikeParameters(SpikeParameters(:,type)==0,invF2)+(rand(sizeUndef,1)-0.5)./20,200,rgb('DimGray'),'.');
%    scatter(SpikeParameters(SpikeParameters(:,type)==1,TtoP)+(rand(sizeExc,1)-0.5)./50,SpikeParameters(SpikeParameters(:,type)==1,invF2)+(rand(sizeExc,1)-0.5)./20,200,rgb('Tomato'),'.');
%    scatter(SpikeParameters(SpikeParameters(:,type)==2,TtoP)+(rand(sizeInh,1)-0.5)./50,SpikeParameters(SpikeParameters(:,type)==2,invF2)+(rand(sizeInh,1)-0.5)./20,200,rgb('RoyalBlue'),'.');
%    xlim([0.1 0.9]);
%    ylim([0.2 1.6]);
%
%  subplot(3,3,8:9)
%    [hi,blips]=hist(SpikeParameters(:,TtoP),[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8]);
%    bar(blips,log(hi));
%  %    hist(SpikeParameters(:,TtoP),[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8]);
%  %    xlim([0.1 0.9]);
%    xlabel('Peak to Trough');
%    ylabel('log(count)');
%  subplot(3,3,[1,4])
%    [h,bins]=hist(SpikeParameters(:,invF2),[0.3:0.1:1.4]);
%    barh(bins,log(h));
%    set(gca,'xdir','r');
%    ylim([0.2 1.6]);
%    ylabel('Spike Width');
%    ylabel('log(count)');
%  struc = input('Structure name?','s')
%  suptitle([struc ' - nCells = ' int2str(size(SpikeParameters,1))])
%
%  %  plot2svg(['SpikeParameters-' struc '.svg'],gcf);
%
%  if ~isempty(clusters)
%    figure;
%    hold on;
%    sizeClu1=size(SpikeParameters(clusters==1,:),1);
%    sizeClu2=size(SpikeParameters(clusters==2,:),1);
%    sizeClu0=size(SpikeParameters(clusters==0,:),1);
%    scatter(SpikeParameters(clusters==0,TtoP)+(rand(sizeClu0,1)-0.5)./50,SpikeParameters(clusters==0,invF2)+(rand(sizeClu0,1)-0.5)./20,200,rgb('DimGrey'),'.');
%    scatter(SpikeParameters(clusters==1,TtoP)+(rand(sizeClu1,1)-0.5)./50,SpikeParameters(clusters==1,invF2)+(rand(sizeClu1,1)-0.5)./20,200,rgb('Tomato'),'.');
%    scatter(SpikeParameters(clusters==2,TtoP)+(rand(sizeClu2,1)-0.5)./50,SpikeParameters(clusters==2,invF2)+(rand(sizeClu2,1)-0.5)./20,200,rgb('RoyalBlue'),'.');
%    xlabel('Peak to Trough');
%    ylabel('Spike Width');
%    xlim([0.1 0.9]);
%    ylim([0.2 1.6]);
%  end
%
%
%  %%% TtoP vs invF2 vs FR (3D)
%  figure;
%  hold on;
%  scatter3(SpikeParameters(SpikeParameters(:,type)==0,TtoP)+rand(sizeUndef,1)./50,SpikeParameters(SpikeParameters(:,type)==0,invF2)+rand(sizeUndef,1)./20,log(SpikeParameters(SpikeParameters(:,type)==0,meanFR)),100,rgb('DimGray'),'.');
%  scatter3(SpikeParameters(SpikeParameters(:,type)==1,TtoP)+rand(sizeExc,1)./50,SpikeParameters(SpikeParameters(:,type)==1,invF2)+rand(sizeExc,1)./20,log(SpikeParameters(SpikeParameters(:,type)==1,meanFR)),300,rgb('Tomato'),'.');
%  scatter3(SpikeParameters(SpikeParameters(:,type)==2,TtoP)+rand(sizeInh,1)./50,SpikeParameters(SpikeParameters(:,type)==2,invF2)+rand(sizeInh,1)./20,log(SpikeParameters(SpikeParameters(:,type)==2,meanFR)),300,rgb('RoyalBlue'),'.');
%  xlabel('Peak to Trough');
%  ylabel('Spike Width');
%  zlabel('Firing Rate');
