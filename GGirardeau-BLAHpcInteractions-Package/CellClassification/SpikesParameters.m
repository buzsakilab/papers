function parameters = SpikesParameters (unitlist, session)

%SpikeParameters - Calculates Spike Parameters needed to determine putative cell types.
%
%  USAGE
%
%    parameters = SpikeParameters (unitlist)
%
%  INPUT
%
%    unitlist           [shank cluster] nx2 matrix
%    session            path to session
%
%  OUTPUT
%
%    parameters		matrix n x [shank unit TtoP invF invF2 RatioTtoP asym meanFR]
%
%  with:
%    TtoP            	Trough to Peak duratio in ms
%    invF	        inverse frequency (=duration) of the spike calculated by Wavelet transform.
%    invF2		inverse frequency (=duration) of the spike calculated by FFT
%    RatioTtoP		Ratio of absolute values peak/trough.
%    asym		Spike aymmetry (Royer 2012)
%    meanFR		Mean firing rate.
%
%  SEE
%
%    See also : GetSpikeWaveforms getWavelet + my_spectrum (Eran Stark) (dependencies)
%
% June 2014, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

cd(session)
xml=session(end-13:end)
SetCurrentSession(xml)

if strcmp(unitlist,'all')
  unitlist=GetUnits;
end


% Initialize variable
parameters=[];

for i = 1: size(unitlist,1)
  unit=(unitlist(i,:));

  %%%%%%%%%%% Trough to Peak duration
  wv=GetSpikeWaveforms(unit);
  allmeanwv=squeeze(mean(wv,1));
  [channel,datapoint]=find(allmeanwv==min(min(allmeanwv)));
  maxmeanwv=allmeanwv(channel,:);
  troughPos=find(maxmeanwv==min(maxmeanwv));
  peakPos=find(maxmeanwv==max(maxmeanwv(troughPos:end)));
  peakPos(peakPos<troughPos)=[];
  if length(peakPos)>1
    peakPos=peakPos(1);
  end
  TtoP=(peakPos-troughPos)*1000/20000;% Peak to trough nBins*1000/20000 -->msec.
  RatioTtoP=abs(min(maxmeanwv))/abs(max(maxmeanwv(troughPos:end)));

  %%%%%%%%%%% Inverse Frequency (duration) using Wavelet transform
  [wave,f,t]=getWavelet(maxmeanwv',20000,500,3000,128,0,0);
  meanWavelet=mean(wave,2);
%    figure;
%    imagesc(t,f,wave);
%    troughWavelet=wave(:,troughPos);
%    set(gca,'YDir','normal');
  SpikeFreq=f(meanWavelet==max(meanWavelet));
  invF=(1/SpikeFreq)*1000;%(= spike width in ms)

  %%%%%%%%%%% Inverse Frequency (duration) using Eran's my_spectrum
  [powspectra,f2]=my_spectrum(maxmeanwv', 1024, 20000, hanning(32), 0, 0 );
  SpikeFreq2=f2(powspectra==max(powspectra));
  invF2=(1/SpikeFreq2)*1000;%(= spike width in ms)

%    %%%%%%%%%% CheckPlot
%    figure;
%    plot(f2,powspectra,'k');
%    hold on;
%    plot(f,meanWavelet,'r');
%    plot(f,troughWavelet,'g')
%    xlim([0 3000]);

  %%%%%%%%%%% Spike asymmetry (Royer 2012)
  firstpeak=max(maxmeanwv(1:troughPos));
  secondpeak=max(maxmeanwv(troughPos:end));
  asym=(firstpeak-secondpeak)/(firstpeak+secondpeak);

  %%%%%%%%%%% Mean FR
  spikes=GetSpikeTimes(unit);
  meanFR=length(spikes)/(spikes(end)-spikes(1));

  %%%%%%%%%% Grouping variable
  parameters=[parameters;[unit] TtoP invF invF2 RatioTtoP asym meanFR];

end

SpikeParameters=parameters;
save([xml '-SpikeParameters'],'SpikeParameters');
cd ..


