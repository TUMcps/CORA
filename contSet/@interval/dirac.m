function res = dirac(varargin)
% dirac - overloaded built-in dirac function
%    note: we do not higher-order derivatives of the dirac function
%
%    elementwise computation:
%       dirac(I)   = [0,0]       if 0 not in I
%                    [0,Inf]     otherwise
%       dirac(1,I) = [0,0]       if 0 not in I
%                    [-Inf,Inf]  if 0 in I
%                    [-Inf,0]    if 0 = min(I)
%                    [0,Inf]     if 0 = max(I)
%
% Syntax:
%    res = dirac(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-2;3],[3;4]);
%    res = dirac(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       17-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 1
    I = varargin{1};
    n = 0;
elseif nargin == 2
    n = varargin{1};
    I = varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% init with zeros
res = interval(zeros(size(I)),zeros(size(I)));

% sign function to check which dimensions contain 0
signI = sign(I);

% set the upper bound of those dimensions to Inf
if n == 0
    res.sup(abs(signI.inf + signI.sup) <= 1) = Inf;
elseif n == 1
    res.inf(signI.inf + signI.sup == 0) = -Inf;
    res.sup(signI.inf + signI.sup == 0) = Inf;
    res.inf(signI.inf + signI.sup == 1) = -Inf;
    res.sup(signI.inf + signI.sup == -1) = Inf;
else
    res.inf(abs(signI.inf + signI.sup) <= 1) = -Inf;
    res.sup(abs(signI.inf + signI.sup) <= 1) = Inf;
end

% ------------------------------ END OF CODE ------------------------------
