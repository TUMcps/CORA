function [result] = in(obj1,obj2)
% in - determines if elements of a zonotope obj2 are in a mptPolytope 
%
% Syntax:  
%    [result] = in(obj1,obj2)
%
% Inputs:
%    obj1 - mptPolytope object
%    obj2 - zonotope object
%
% Outputs:
%    result - 1/0 if zonotope is in, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      12-February-2012
% Last update:  27-July-2016
%               03-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%obtain the C matrix and d vector
try %MPT3
    C = obj1.P.A;
    d = obj1.P.b;
catch %MPT 2
    [C,d] = double(obj1.P); 
end

%affine map
Znew = C*obj2 + (-d);

%compute interval
IHnew = interval(Znew);

if all(infimum(IHnew)<=0)
    result = 1;
else
    result = 0;
end


%------------- END OF CODE --------------
