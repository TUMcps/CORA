function [value] = sup(Z)
% sup - Determines $sup(||x||_\infty),x in Z$, whereas sup is the operator
% determining the supremum of its argument.
%
% Syntax:  
%    [value]=sup(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    value - supremum of $||x||_\infty,x in Z$
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    supremum=sup(Z)-->supremum=??
%
% Other m-files required: intervalhull
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
I=interval(Z);

%determine vector with greatest infinity norm within the interval hull
N1=norm(infimum(I),inf);
N2=norm(supremum(I),inf);
value=max(N1,N2);

%------------- END OF CODE --------------