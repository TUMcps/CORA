function I = enlarge(I,factor)
% enlarge - Enlarges an interval object around its center
%
% Syntax:
%    obj = enlarge(I,factor)
%
% Inputs:
%    I - interval object
%    factor - enlarging factor (scalar or column vector)
%
% Outputs:
%    I - enlarged interval object
%
% Example: 
%    I = interval([-1;-2],[3;4]);
%    I = enlarge(I,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       22-July-2016 
% Last update:   28-August-2019
%                03-December-2023 (MW, support unbounded sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get center and radius
c = center(I);
r = rad(I);

% enlarged intervals
I.inf = c - r.*factor; 
I.sup = c + r.*factor;

% unbounded dimensions: r is Inf
if any(r == Inf & factor < 1)
    % if factor < 1, there is no center -> return error (NaN)
    throw(CORAerror('CORA:notSupported'));
else
    % factor > 1, dimension expands to [-Inf,Inf]
    Infdims = r == Inf & factor > 1;
    I.inf(Infdims) = -Inf;
    I.sup(Infdims) = Inf;
end

% ------------------------------ END OF CODE ------------------------------
