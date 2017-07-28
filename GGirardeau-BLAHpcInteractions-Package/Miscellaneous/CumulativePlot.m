function CumulativePlot(values,varargin)

%Cumulative Plot - Draws a cumulative plot for a distribution
%
%  USAGE
%
%    FunctionName (x,<options>)
%
%    values		list of values (vector)                   
%    <options>      	optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties            Values
%    -------------------------------------------------------------------------
%     'normalize'	'on','off' (default)
%     'newfig'		'on' (default),'off'
%     'color'		nem of color as in rgb 'BrightOrange'
%    =========================================================================
%
%  OUTPUT
%
%    figure    
%
%  NOTE
%     NaN values are discarded
%
%  SEE
%
%    See also 
%
% Gabrielle Girardeau, January 2016
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults
newfig = 'on';
normalize = 'off';
color = 'Black';

% Check number of inputs
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).');
end

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
		case 'newfig',
			newfig = varargin{i+1};
		case 'normalize',
			normalize = lower(varargin{i+1});
		case 'color'
			color = varargin{i+1};
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FunctionName">FunctionName</a>'' for details).']);
	end
end

OrdValues=sort(values);
OrdValues(isnan(OrdValues))=[];
summ=1:length(OrdValues);

if strcmp(normalize,'on')
  summ=summ./length(OrdValues);
end

if strcmp(newfig,'on')
  figure;
end
plot(OrdValues,summ,'Color',rgb(color),'LineWidth',2);




