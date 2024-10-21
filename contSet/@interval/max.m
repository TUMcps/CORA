function res = max(I, Y, varargin)
% max - computes the maximum of I and Y
%
% Syntax:
%    res = max(I, Y)
%
% Inputs:
%    I - interval object
%    Y - interval or numeric 
%    varargin - additional parameters for build-in max function
%
% Outputs:
%    res - interval
%
% Example: 
%    I = interval([-2;-1],[2;1]);
%    res = max(I, 0)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: max, interval/supremum

% Authors:       Tobias Ladner
% Written:       16-December-2022
% Last update:   11-April-2024 (TL, single input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    % return supremum
    res = I.sup;
    return
end

% parse input
inputArgsCheck({ ...
    {I, 'att', 'interval'};
    {Y, 'att', {'contSet', 'numeric'}};
});

if isnumeric(Y)
    % check dimensions
    if ~isempty(Y) && ~isscalar(Y) && ~all(dim(I) == size(Y), "all")
        throw(CORAerror('CORA:dimensionMismatch', I, Y));
    end
    res = interval(max(I.inf, Y, varargin{:}), max(I.sup, Y, varargin{:}));
    return;
end

if ~isa(Y, 'interval')
    % convert contSet to interval
    Y = interval(Y);
end

res = interval( ...
    max(I.inf, Y.inf, varargin{:}), ...
    max(I.sup, Y.sup, varargin{:}) ...
);
    
   
% ------------------------------ END OF CODE ------------------------------
