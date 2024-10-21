function [infi,isMin] = infimum(obj)
% infimum - returns the infimum of an STL interval
%
% Syntax:
%    infi = infimum(obj)
%    [infi,isMin] = infimum(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    infi - infimum of interval
%    isMin - boolean indicating if infimum is included in interval
%
% Example: 
%    I = stlInterval(1,2,true,true);
%    [infi,isMin] = infimum(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

infi = obj.lower;
isMin = obj.leftClosed;

% ------------------------------ END OF CODE ------------------------------
