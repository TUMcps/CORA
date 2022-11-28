function Z = convHull(Z,S)
% convHull - computes an enclosure for the convex hull of a zonotope and
%    another set or a point
%
% Syntax:  
%    Z = convHull(Z,S)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%
% Outputs:
%    Z - zonotope enclosing the convex hull
%
% Example: 
%    Z1 = zonotope([2 1 0; 2 0 1]);
%    Z2 = zonotope([-2 1 0; -2 0 1]);
%
%    Z = convHull(Z1,Z2);
%
%    figure; hold on;
%    plot(Z1,[1,2],'FaceColor','r');
%    plot(Z2,[1,2],'FaceColor','b');
%    plot(Z,[1,2],'g','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:        Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 1
    return;
end 
if isempty(S)
    return;
end
% determine zonotope object
if ~isa(Z,'zonotope')
    temp = Z;
    Z = S;
    S = temp;
end

% different cases depending on the class of the second summand
if isa(S,'zonotope') || isa(S,'interval') || isnumeric(S)

    Z = enclose(Z,zonotope(S));

elseif isa(S,'mptPolytope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'polyZonotope') || ...
       isa(S,'conPolyZono')

    Z = convHull(S,Z);        

else

    % throw error for given arguments
    throw(CORAerror('CORA:noops',Z,S));

end

%------------- END OF CODE --------------