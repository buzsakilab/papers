function  allparam=SpikeParametersAll(ratnumber)

%SpikeParametersAll - Merges SpikeParameters matrices across sessions for each rat. 
%
%  USAGE
%
%    allparams = SpikeParametersAll (ratnumber)
%
%  INPUT
%
%    ratnumber          rat to process
%
%  OUTPUT
%
%    allparams		matrix n x [rat sess shank unit shank unit TtoP invF invF2 RatioTtoP asym meanFR ccgclassif]
%
%  with:
%    rat
%    session
%    shank
%    unit
%    TtoP            	Trough to Peak duratio in ms
%    invF	        inverse frequency (=duration) of the spike calculated by Wavelet transform.
%    invF2		inverse frequency (=duration) of the spike calculated by FFT
%    RatioTtoP		Ratio of absolute values peak/trough.
%    asym		Spike aymmetry (Royer 2012)
%    meanFR		Mean firing rate.
%    MonoSynType        Cell type as established by the monosynaptic connections (CCGs)
%
%  SEE
%
%    See also : SpikeParameters
%
% June 2014, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/sessionindexing.mat')

allparam=[];
sessionlist=xmlpath(ratsessionindex(:,1)==ratnumber);
sessionnumbers=ratsessionindex(ratsessionindex(:,1)==ratnumber,2);
for i=1:size(sessionlist)
  currentsession=sessionlist{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  load([xml '-SpikeParameters.mat']);
  load([xml '-MonoSynConvClick.mat']);
  type=zeros(size(SpikeParameters,1),1);
  if ~isempty(FinalExcCellList)
    type(ismember(SpikeParameters(:,1:2),FinalExcCellList,'rows'))=1;
  end
  if ~isempty (FinalInhCellList)
    type(ismember(SpikeParameters(:,1:2),FinalInhCellList,'rows'))=2;
  end
  sess=ones(size(SpikeParameters,1),1).*sessionnumbers(i);
  sessionparam=[sess SpikeParameters type];
  allparam=[allparam;sessionparam];
end



