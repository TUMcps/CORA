function msg = intervalNaNInf()
% intervalNaNInf - standardized format if a mathematical function
%    of @interval returns NaN and/or Inf
%
% Syntax:  
%    msg = intervalNaNInf()
%
% Inputs:
%    -
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
filename = st(2).name; % calling function name

msg = sprintf(...
    "The interval is not inside the valid domain of the function %s().",...
    filename);

%------------- END OF CODE --------------
