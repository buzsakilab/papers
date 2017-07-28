function CreateCellIndex(session)
% Creates an index of all units for a single session

load('/media/Data-01/All-Rats/sessionindexing.mat');

cd(session);
xml=session(end-13:end);
SetCurrentSession(xml);

%%%%%%%%%%%%%%%%%%%%%%% Load Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get All Spikes
spikes1=GetSpikeTimes('output','numbered');
spikes2=GetSpikeTimes('output','full');
spikes=[spikes2 spikes1(:,2)];
spikes(spikes(:,3)==0,:)=[];
spikes(spikes(:,3)==1,:)=[];

% Renumber IDs after removing MUA and artifacts (spikes : t / shank / cell / ID)
ids=unique(spikes(:,4));
for ii=1:length(ids)
  spikes(spikes(:,4)==ids(ii),4)=ii;
end
if sum(diff(unique(spikes(:,4)))-1)~=0
  warning('Problem with IDs');
end
Index=unique(spikes(:,2:4),'rows'); % Build index for session

ratsess=ratsessionindex(strcmp([session '/'],xmlpath),1:2);

Index=[repmat(ratsess,length(Index),1) Index];
save([xml '-CellIndex.mat'],'Index');