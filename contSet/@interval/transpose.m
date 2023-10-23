function I = transpose(I)
% transpose - Overloaded '.'' operator for single operand
%
% Syntax:
%    I = transpose(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([-2;-3;-4],[5;6;7]);
%    I.'
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       07-February-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I.inf = I.inf.';
I.sup = I.sup.';

% ------------------------------ END OF CODE ------------------------------
