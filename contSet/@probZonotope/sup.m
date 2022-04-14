function value = sup(probZ)
% sup - Determines $sup(||x||_\infty),x in Z$, whereas sup is the operator
%    determining the supremum of its argument.
%
% Syntax:  
%    value = sup(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    value - supremum of $||x||_\infty,x in Z$
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    sup(probZ)
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      30-September-2006 
% Last update:  22-March-2007
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%convert zonotope to interval
Int=interval(probZ);

%determine vector with greatest infinity norm within the interval hull
N1=norm(infimum(Int),Inf);
N2=norm(supremum(Int),Inf);
value=max([N1,N2]);

%------------- END OF CODE --------------