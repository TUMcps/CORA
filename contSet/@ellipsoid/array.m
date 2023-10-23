function E = array(varargin)
% array - returns an array of empty ellipsoids of specified size
%
% Syntax:
%    c = array(n)
%    c = array(n1,n2,...)
%
% Inputs:
%    n   - number of rows and columns
%    ni  - length of i-th dimension
%
% Outputs:
%    E - ellipsoid array of specified size
%
% Example: 
%    E = ellipsoid.array(2,3);
%    size(E) % ans = [2,3]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

E = emptyClassArray('ellipsoid',varargin{:});

% ------------------------------ END OF CODE ------------------------------
