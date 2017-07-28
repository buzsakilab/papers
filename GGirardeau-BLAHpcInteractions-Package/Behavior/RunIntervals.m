function runs=RunIntervals(index)

% index = vector of subsession number that you want the intervals from

beg=GetEvents({'beg.*'});
endd=GetEvents({'end.*'});
evts=[beg endd];

runs=[];

for i=1:size(index,2)
  runs=[runs;evts(index(i),:)];
end


