function I = ctranspose(I)
% ctranspose - Overloaded ''' operator for single operand
%
% Syntax:
%    I = ctranspose(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([-1;-2],[3;4]);
%    I'
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       14-August-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I.inf = I.inf.';
I.sup = I.sup.';

% ------------------------------ END OF CODE ------------------------------
