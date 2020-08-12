function [E] = plus(E1,E2,L)
% plus - Overloaded '+' operator for the addition of two ellipsoids or
% ellipsoid and vector
%
% Syntax:  
%    [E] = plus(E1,E2)
%
% Inputs:
%    E1          - Ellipsoid object/vector
%    E2          - Ellipsoid object/vector
%    (optional)L - Unit directions over which length(N) different
%    lplus(E1,E2,L(:,i)) results are intersected
%
% Outputs:
%    E - Ellipsoid after addition of two ellipsoids/ellipsoid and vector
%
% Example: 
%    E1=ellipsoid([1 0; 0 1]);
%    E2=ellipsoid([1 1; 1 1]);
%    E =E1+E2;
%    n = length(E1.Q);
%    X = randn(n,N);
%    L = X./repmat(sqrt(sum(X.^2,1)),n,1);
%    E_alternative = plus(E1,E2,L);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if (~isa(E1,'ellipsoid') && ~isa(E1,'double')) || ...
   (~isa(E1,'ellipsoid') && ~isa(E1,'double'))
    error('Wrong type of arguments');
end
if isa(E1,'double') && isa(E2,'double')
    error('At least one argument has to be of class "ellipsoid"');
end
if isa(E1,'double')
    if size(E1,2)~=1 || (size(E1,1)~=1 && size(E1,1)~=size(E2.q,1))
        error('Argument has to be either a scalar or a vector with appropriate length');
    end
    E = ellipsoid(E2.Q,E2.q+E1);
    return;
elseif isa(E2,'double')
    if size(E2,2)~=1 || (size(E2,1)~=1 && size(E2,1)~=size(E1.q,1))
        error('Argument has to be either a scalar or a vector with appropriate length');
    end
    E = ellipsoid(E1.Q,E1.q+E2);
    return;
end
if E1.isdegenerate || E2.isdegenerate
    error('Ellipsoids have to be full-dimensional');
end
n = E1.dim;
%compute N random directions
if ~exist('L','var')
    %idea is that the intersection has smaller volume; if not, use 
    %use equal partition; compute enough directions such that each
    %dimension is "covered"
    L = eq_point_set(n-1,2*n);
    counter = 1;
    while rank(L)<n
        L = eq_point_set(n-1,2*n+counter);
        counter = counter + 1;
    end
end
E = lplus(E1,E2,L(:,1));
for i=1:length(L)
    Etmp = lplus(E1,E2,L(:,i));
    %intersection is necessarily nonempty
    E = and(E,Etmp,false);
end
%------------- END OF CODE --------------