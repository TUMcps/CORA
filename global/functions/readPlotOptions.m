function [linespec,NVpairs] = readPlotOptions(plotOptions)
% readPlotOptions - separates LineSpecification options from
%    name-value pairs in plot options
%
% Syntax:  
%    [linespec,NVpairs] = readPlotOptions(plotOptions)
%
% Inputs:
%    plotOptions - LineSpecification options + Name-Value pairs
%
% Outputs:
%    linespec - LineSpecification options for Matlab plots
%    NVpairs  - Name-Value pairs for Matlab plots
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/plot (all classes)

% Author:        Mark Wetzlinger
% Written:       14-July-2020 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% default values
linespec = 'b';

% process linespec and Name-Value pairs
allowedChars = '-:.+o*xsd^v><phrgbcmykw';

% determine start of name-value pairs
if all(ismember(plotOptions{1},allowedChars))
    % linespec is given, therefore name-value pairs start afterwards
    linespec = plotOptions{1}; startNV = 2;
else
    % no linespec given, only name-value pairs
    startNV = 1;
end

% remaining inputs are name-value pairs
NVpairs = plotOptions(startNV:end);

end

%------------- END OF CODE --------------
