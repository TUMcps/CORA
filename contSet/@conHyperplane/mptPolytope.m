function P = mptPolytope(obj)
% mptPolytope - Converts a constrained hyperplane to a mptPolytope object
%
% Syntax:  
%    P = mptPolytope(obj)
%
% Inputs:
%    obj - conHyperplane object
%
% Outputs:
%    P - mptPolytope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/mptPolytope

% Author:       Niklas Kochdumper
% Written:      26-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(obj.C)
        A = [obj.h.c';-obj.h.c'];
        b = [obj.h.d;-obj.h.d];       
    else
        A = [obj.h.c';-obj.h.c';obj.C];
        b = [obj.h.d;-obj.h.d;obj.d];
    end

    P = mptPolytope(A,b);

end

%------------- END OF CODE --------------