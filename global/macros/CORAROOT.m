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

% Authors:       ???, Tobias Ladner
% Written:       ---
% Last update:   23-May-2025 (TL, major speedup)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

s = mfilename("fullpath");
corapath = fileparts(fileparts(fileparts(s)));

% ------------------------------ END OF CODE ------------------------------
