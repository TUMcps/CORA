function P = mptPolytope(obj)
% mptPolytope - Converts an halfspace object to a polytope object
%
% Syntax:  
%    P = mptPolytope(obj)
%
% Inputs:
%    obj - halfspace object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    hs = halfspace(C,d);
%    P = mptPolytope(hs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/mptPolytope

% Author:       Victor Charlent
% Written:      28/06/2016
% Last update:  17-March-2017, Matthias Althoff
% Last revision:---

%------------- BEGIN CODE --------------

    A = obj.c';
    b = obj.d;

    P = mptPolytope(A,b);

end

%------------- END OF CODE --------------

