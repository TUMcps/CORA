function corapath = CORAROOT()
% CORAROOT - returns the CORA root path
%
% Syntax:
%    corapath = CORAROOT()
%
% Inputs:
%    -
%
% Outputs:
%    corapath (string) - directory of CORA
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

s = which('coraroot');
corapath = fileparts(fileparts(fileparts(s)));

% ------------------------------ END OF CODE ------------------------------
