function varargout = LinParamWrapper_initReach(varargin)
% LinParamWrapper_initReach - wrapper for initReach function from linParamSys
%
% Syntax:
%    varargout = LinParamWrapper_initReach(varargin)
%
% Inputs:
%    -
%
% Outputs:
%    ?
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[varargout{1:nargout}] = initReach(varargin{:});

% ------------------------------ END OF CODE ------------------------------
