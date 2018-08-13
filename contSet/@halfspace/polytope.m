function [P]=polytope(obj)
% polytope - Converts an halfspace object to a polytope object
%
% Syntax:  
%    [P]=polytope(obj)
%
% Inputs:
%    obj - halfspace object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    hs=halfspace(C,d);
%    P=polytope(hs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor Charlent
% Written:      28/06/2016
% Last update:  17-March-2017, Matthias Althoff
% Last revision:---

%------------- BEGIN CODE --------------

A = get(obj, 'c')';
b = get(obj, 'd');

P = mptPolytope(A,b);

end


%------------- END OF CODE --------------

