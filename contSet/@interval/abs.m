function res = abs(I)
% abs - returns the absolute value of an interval
%
% Syntax:
%    res = abs(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-2;-1],[3;4]);
%    res = abs(I);
%
%    figure; hold on;
%    plot(I,[1,2],'b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2015
% Last update:   14-February-2015
%                12-October-2015
%                27-May-2024 (MW, simplify and speed up code)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init resulting interval object without going through constructor (fast)
res = I;

% use size(I.inf) since size(I) is a bit slow... include sparse handling
res.inf = max(cat(3, zeros(size(I.inf)), full(I.inf), full(-I.sup)), [], 3);
res.sup = max(cat(3, full(-I.inf), full(I.sup)), [], 3);

% ------------------------------ END OF CODE ------------------------------
