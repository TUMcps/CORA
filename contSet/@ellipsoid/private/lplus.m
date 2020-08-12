function [E] = lplus(E1,E2,l)
% lplus - Computes the Minkowski sum of two ellipsoids such that resulting
% overapproximation given by E is tight in direction l
%
% Syntax:
%    [E] = lplus(E1,E2,l)
%
% Inputs:
%    E1 - Ellipsoid object
%    E2 - Ellipsoid object 
%    l  - unit direction
%
% Outputs:
%    E - Ellipsoid after addition of two ellipsoids
%
% Example: 
%    E1=ellipsoid([1 0; 0 1]);
%    E2=ellipsoid([1 1; 1 1]);
%    l =[1;0];
%    E =lplus(E1,E2,l);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus
%
% References:
%    [1] https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
q = E1.q+E2.q;
Q1 = E1.Q;
Q2 = E2.Q;
%for details, see [1]
Q = (sqrt(l'*Q1*l)+sqrt(l'*Q2*l))*(Q1/sqrt(l'*Q1*l)+Q2/sqrt(l'*Q2*l));
E = ellipsoid(Q,q);
%------------- END OF CODE --------------