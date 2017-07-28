function pos = GetCleanPositions(varargin)
%GetCleanPositions - Gets positions (1 LED - discards blue LED) in real coordinates with pixel = 0.43cm (for linear track) 

% Defaults
show = 'on' ;

% Parse options
for i = 1:2:length(varargin),
	switch(lower(varargin{i})),
		case 'show',
			show = varargin{i+1};
	end
end

pos=GetPositions('coordinates','real','pixel',0.43,'discard','none');
CleanPos(pos);
pos(:,2:3)=[];
pos=InterpolateNaNPos(pos);
if strcmp(show,'on')
  figure;
  PlotXY(pos(:,1:2));
end