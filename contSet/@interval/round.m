function I = round(I,varargin)
% round - rounds each element of the interval to the given precision
%
% Syntax:
%    I = round(I)
%    I = round(I,N)
%
% Inputs:
%    I - interval object
%    N - number of digits
%
% Outputs:
%    I - rounded interval object
%
% Example: 
%    I = interval([-1.5;1.1],[1.2;2.4]);
%    I_rounded = round(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: round

% Authors:       Tobias Ladner
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
N = setDefaultValues({0},varargin);
inputArgsCheck({ ...
    {I,'att','interval'}, ...
    {N,'att','numeric',{'scalar','nonnegative','finite'}} ...
})

% round
I.inf = round(I.inf,N,varargin{2:end});
I.sup = round(I.sup,N,varargin{2:end});

end

% ------------------------------ END OF CODE ------------------------------
