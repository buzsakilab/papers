function [allccgallrats,Idx] = AllRippleCCGsAll

%AllRippleCCGs_2 - Merges all spikes-ripples CCGs for across animals and sessions.
%
%  USAGE
%
%    [allccgallrats,Idx] = AllRippleCCGsAll
%
%  OUTPUT
%
%    allccgallrats           matrix of cross-correolgrams between spikes and ripples (1 line/cell) for all cells across animals and sessions
%    Idx                     Cell Index
%
%  NOTE
%
%  SEE
%
%    See also : AllRippleCCGs_2
%
% Gabrielle Girardeau, 2016
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

load('/media/Data-01/All-Rats/sessionindexing.mat')

allccgallrats.ccg1=[];
allccgallrats.ccg2=[];
allccgallrats.gain1=[];
allccgallrats.gain2=[];

allccgallrats.pre.ccg1=[];
allccgallrats.pre.ccg2=[];
allccgallrats.pre.gain1=[];
allccgallrats.pre.gain2=[];

allccgallrats.post.ccg1=[];
allccgallrats.post.ccg2=[]
allccgallrats.post.gain1=[];
allccgallrats.post.gain2=[];

Idx.all=[];
Idx.prepost=[];

for i=1:size(xmlpath)
  currentsession=xmlpath{i}
  cd(currentsession);
  xml=currentsession(end-14:end-1);
  if exist([xml '-RippleCCGs.mat'])==2
    load([xml '-RippleCCGs.mat']);
    allccgallrats.ccg1=[allccgallrats.ccg1;allccg.ccg1];
    allccgallrats.ccg2=[allccgallrats.ccg2;allccg.ccg2];
    allccgallrats.gain1=[allccgallrats.gain1;allccg.gain1];
    allccgallrats.gain2=[allccgallrats.gain2;allccg.gain2];
    ratsess=repmat(ratsessionindex(i,:),[length(idx) 1]);
    Idx.all=[Idx.all;[ratsess idx]];
    if isfield(allccg,'pre') 
      allccgallrats.pre.ccg1=[allccgallrats.pre.ccg1;allccg.pre.ccg1];
      allccgallrats.pre.ccg2=[allccgallrats.pre.ccg2;allccg.pre.ccg2];
      allccgallrats.pre.gain1=[allccgallrats.pre.gain1;allccg.pre.gain1];
      allccgallrats.pre.gain2=[allccgallrats.pre.gain2;allccg.pre.gain2];
      allccgallrats.post.ccg1=[allccgallrats.post.ccg1;allccg.post.ccg1];
      allccgallrats.post.ccg2=[allccgallrats.post.ccg2;allccg.post.ccg2];
      allccgallrats.post.gain1=[allccgallrats.post.gain1;allccg.post.gain1];
      allccgallrats.post.gain2=[allccgallrats.post.gain2;allccg.post.gain2];
      Idx.prepost=[Idx.prepost;[ratsess idx]];
      allccgallrats.t1=allccg.t1;
      allccgallrats.t2=allccg.t2;
    end
  end
end

binSize1=0.001;
binSize2=0.01;
duration1=0.4;
duration2=4;
idxall=Idx;
allccg=allccgallrats;
cd('/media/Data-01/All-Rats/');

save('AllRats-AllRippleCCGs','idxall','allccg','binSize1','binSize2','duration1','duration2');