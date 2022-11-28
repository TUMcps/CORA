function res = or(I,S)
% or - computes the union of an interval and a set
%
% Syntax:  
%    res = or(I,S)
%
% Inputs:
%    I - interval object
%    S - contSet object
%
% Outputs:
%    res - union of the two sets
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

% Author:       Niklas Kochdumper
% Written:      25-July-2019
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------
if isempty(S) && ~isempty(I)
    res = I; return;
end
    
% determine the interval object
if ~isa(I,'interval')
    temp = I;
    I = S;
    S = temp;
end

% different cases depending on the class of the summand
if isa(S,'interval')
    
    if dim(I) ~= dim(S)
        throw(CORAerror('CORA:dimensionMismatch',I,S));
    end

    res = interval(min([I.inf,S.inf],[],2), ...
                   max([I.sup,S.sup],[],2));

elseif isnumeric(S)

    res = interval(min([I.inf,S],[],2), ...
                   max([I.sup,S],[],2));

elseif isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'polyZonotope') || ...
       isa(S,'mptPolytope') || isa(S,'conPolyZono')

    res = S | I;

else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',I,S));
end

%------------- END OF CODE --------------