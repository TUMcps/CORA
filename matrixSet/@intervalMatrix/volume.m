function vol = volume(matI)
% volume - computes the volume of an interval matrix by computing the 
%    volume of the corresponding interval
%
% Syntax:  
%    vol = volume(matI)
%
% Inputs:
%    matI - interval matrix
%
% Outputs:
%    vol - volume
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      24-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to a interval
I = interval(matI);

%compute volume of the interval
vol = volume(I);

%------------- END OF CODE --------------