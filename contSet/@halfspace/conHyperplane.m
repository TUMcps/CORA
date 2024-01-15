function hyp = conHyperplane(hs,varargin)
% conHyperplane - conversion of halfspace to conHyperplane
%
% Syntax:
%    hyp = conHyperplane(hs)
%    hyp = conHyperplane(hs,C,d)
%
% Inputs:
%    hs - halfspace object
%    C - constraint matrix
%    d - constraint vector
%
% Outputs:
%    hyp - conHyperplane object
%
% Example:
%    a = [-0.5; -1; 0.1]; b = -1;
%    hs = halfspace(a,b);
%    C = [-0.6 0.8 -1.7;...
%          0.6 0.5 -0.8];
%    d = [1; 0.5];
%    hyp = conHyperplane(a,b,C,d);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/conHyperplane

% Authors:       Mark Wetzlinger
% Written:       23-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

hyp = conHyperplane(hs.c,hs.d,varargin{:});

% ------------------------------ END OF CODE ------------------------------
