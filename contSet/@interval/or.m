function S_out = or(I,S,varargin)
% or - computes the union of an interval and a set or point
%
% Syntax:
%    S_out = I | S
%    S_out = or(I,S)
%    S_out = or(I,S,mode)
%
% Inputs:
%    I - interval object
%    S - contSet object
%    mode - 'exact', 'outer', 'inner'
%
% Outputs:
%    S_out - union of the two sets
%
% Example: 
%    I1 = interval([-2;-2],[-1;-1]);
%    I2 = interval([0;0],[2;2]);
%    res = I1 | I2;
%
%    figure; hold on; xlim([-3,3]); ylim([-3,3]);
%    plot(I1,[1,2],'FaceColor','g');
%    plot(I2,[1,2],'FaceColor','b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/split

% Authors:       Niklas Kochdumper
% Written:       25-July-2019
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
mode = setDefaultValues({'outer'},varargin);

% ensure that numeric is second input argument
[I,S] = reorderNumeric(I,S);

% check dimensions
equalDimCheck(I,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < I.precedence
    S_out = or(S,I,varargin{:});
    return
end

% empty set case: union is the other set
if representsa_(S,'emptySet',eps)
    S_out = I;
    return
elseif representsa_(I,'emptySet',eps)
    S_out = S;
    return
end

% interval-interval case
if isa(S,'interval')
    S_out = interval(min(I.inf,S.inf),max(I.sup,S.sup));
    return
end

% numeric
if isnumeric(S)
    S_out = interval(min(I.inf,S),max(I.sup,S));
    return
end

throw(CORAerror('CORA:noops',I,S));

% ------------------------------ END OF CODE ------------------------------
