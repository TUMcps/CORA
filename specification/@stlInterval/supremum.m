function [sup,isMax] = supremum(obj)
% supremum - returns the supremum of an STL interval
%
% Syntax:
%    sup = supremum(obj)
%    [sup,isMax] = supremum(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    sup - supremum of interval
%    isMax - boolean indicating if supremum is included in interval
%
% Example: 
%    I = stlInterval(1,2,true,true);
%    [sup,isMax] = supremum(I)
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

sup = obj.upper;
isMax = obj.rightClosed;

% ------------------------------ END OF CODE ------------------------------
