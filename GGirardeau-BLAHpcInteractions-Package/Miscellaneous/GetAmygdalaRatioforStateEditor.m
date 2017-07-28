function ratio = GetAmygdalaRatioforStateEditor(lfpchannel)

% Calculates and resamples at 1Hz the gamma/low ratio for sleep scoring with amygdala signal.
% To be loaded in The State Editor.

lfp=GetLFP(lfpchannel);
[spectro,t,f]=MTSpectrogram(lfp,'show','off','range',[0 100],'window',1);
bands=SpectrogramBands(spectro,f);
rawratio=[t bands.ratios.amygdala];
timestamps=0:1:round(t(end));
[interpolated,discarded]=Interpolate(rawratio,timestamps);
ratio=interpolated';
save('AmygdalaRatio.mat','ratio');
%  PlotXY(interpolated);
