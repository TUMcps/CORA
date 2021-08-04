function msg = noops(varargin)
% noops - standardized format if a set operation is not implemented
%
% Syntax:  
%    noops(varargin)
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
% See also: 

% Author:        Mark Wetzlinger
% Written:       05-May-2020 
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
    "The function %s is not implemented for the following arguments:\n  %s",...
    filename,classlist);

%------------- END OF CODE --------------
