function S_out = convHull_(I,S,varargin)
% convHull_ - computes the convex hull of an interval and another set or a
%    point
%
% Syntax:
%    res = convHull_(I,S)
%
% Inputs:
%    I - interval object
%    S - contSet object
%
% Outputs:
%    S_out - convex hull of I and S
%
% Example:
%    I = interval([-1;-2],[3;4]);
%    Z = zonotope(ones(2,1),[1 0.5; -0.2 -2]);
%    res = convHull(I,Z);
%
%    figure; hold on;
%    plot(I,[1,2],'b');
%    plot(Z,[1,2],'b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
%                12-March-2021 (MW, add empty case)
% Last revision: 29-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% intervals are convex
if nargin == 1
    S_out = I;
    return;
end

% ensure that numeric is second input argument
[I,S] = reorderNumeric(I,S);

% check dimensions
equalDimCheck(I,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < I.precedence
    S_out = convHull(S,I,varargin{:});
    return
end

% empty set cases
if representsa_(I,'emptySet',eps)
    S_out = S;
    return
elseif representsa_(S,'emptySet',eps)
    S_out = I;
    return
end

% call or function
S_out = I | S;

% ------------------------------ END OF CODE ------------------------------
