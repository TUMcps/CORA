function p = randPoint_(O,N,type,varargin)
% randPoint_ - generates random points within an empty set
%
% Syntax:
%    p = randPoint_(O)
%    p = randPoint_(O,N)
%    p = randPoint_(O,N,type)
%    p = randPoint_(O,'all','extreme')
%
% Inputs:
%    O - emptySet object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme')
%
% Outputs:
%    p - random point in empty set
%
% Example: 
%    O = emptySet(2);
%    p = randPoint(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% unfortunately, we cannot concatentate empty column vectors...
p = double.empty(O.dimension,0);

% ------------------------------ END OF CODE ------------------------------
