function idx = OrderedImagesc(matrix,varargin)

%OrderedImagsc - Plots an Imagesc of a matrix with lines ordered by value.
%
%  USAGE
%
%    idx = OrderedImagesc (matrix,<options>)
%
%    matrix             martix, 1 line per cell/list of values              
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'x'		xvalues (ex:bins)
%     'y'		yvalues
%     'ordertype'     	Order by position of max ('maxi' : default) or minimum ('mini') value or not ('none')
%     'norm'		Normalize values per line (between [0 1]). 'on'(default) or 'off'
%     'idx'		Provide index for ordering (vector same length as matrix) default [];
%     'newfig'		Plot in new figure 'on' (default) 'off'
%    =========================================================================
%
%  OUTPUT
% 
%    figure  
%    idx	reordering index
%
%  NOTE
%
%  
%  SEE
%
%    See also : imagesc.m, ZeroToOne.m
%
% Gabrielle Girardeau Nov 2015
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
x = [];
y = [] ;
norm = 'on';
ordertype = 'maxi';
idx=[];
newfig='on';

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
		case 'x',
			x = varargin{i+1};
		case 'y',
			y = varargin{i+1};
		case 'norm',
			norm = varargin{i+1};
		case 'ordertype',
			ordertype = varargin{i+1};
		case 'idx'
			idx = varargin{i+1};
		case 'newfig'
			newfig = varargin{i+1};	
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end


if strcmp(norm,'on')
  for i = 1: size(matrix,1)
    normatrix(i,:)=ZeroToOne(matrix(i,:));
  end
else
  normatrix=matrix;
end

if ~isempty(idx)
  orderedMatrix=normatrix(idx,:);
else
  if strcmp(ordertype,'maxi')
    % Reorder :position of maximum value
    a=(max(normatrix,[],2));
    for i=1:size(normatrix,1)
      maxpos1(i)=find(normatrix(i,:)==a(i),1);
    end
    [sorted,idx]=sort(maxpos1);
    orderedMatrix=normatrix(idx,:);
  elseif strcmp(ordertype,'mini')
    % Reorder :position of minimum value
    a=(min(normatrix,[],2));
    for i=1:size(normatrix,1)
      minpos(i)=find(normatrix(i,:)==a(i),1);
    end
    [sorted,idx]=sort(minpos);
    orderedMatrix=normatrix(idx,:);
  elseif strcmp(ordertype,'none')
    orderedMatrix=normatrix;
    idx=1:size(normatrix,1);
  end
end

if strcmp(newfig,'on')
  figure;
end
imagesc(x,y,orderedMatrix);



