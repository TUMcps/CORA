function c = center(P)
% center - Returns the chebychev center of a polytope
%
% Syntax:  
%    c = center(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    c - center of the polytope
%
% Example:
%    P = mptPolytope.generateRandom(2);
%    c = center(P);
%
%    figure; hold on;
%    plot(P);
%    plot(c(1),c(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      30-October-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    temp = chebyCenter(P.P);
    c = temp.x;

%------------- END OF CODE --------------