function msg = errNoExactAlg(varargin)
% errNoExactAlg - standardized format if no exact algorithm is implemented
%
% Syntax:  
%    msg = errNoExactAlg(varargin)
%
% Inputs:
%    varargin - sets on which operation was to be performed
%
% Outputs:
%    msg - error message
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: noops

% Author:        Niklas Kochdumper
% Written:       05-February-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% name of operation
st = dbstack;
filename = st(2).file; % calling function name

% read out classes
classlist = "";
for i=1:length(varargin)-1
    classlist = classlist + string(class(varargin{i})) + ", ";
end
% last entry
classlist = classlist + string(class(varargin{end})) + ".";

msg = sprintf(...
    "For the function %s no exact algorithm is implemented for the following arguments:\n  %s",...
    filename,classlist);

%------------- END OF CODE --------------