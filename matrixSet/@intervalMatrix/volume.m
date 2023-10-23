function vol = volume(intMat)
% volume - computes the volume of an interval matrix by computing the 
%    volume of the corresponding interval
%
% Syntax:
%    vol = volume(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
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

% Authors:       Matthias Althoff
% Written:       24-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%conversion to a interval
I = interval(intMat);

%compute volume of the interval
vol = volume(I);

% ------------------------------ END OF CODE ------------------------------
