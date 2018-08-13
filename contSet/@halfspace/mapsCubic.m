function [k,L,Q,C_c,C_g] = maps(obj, A, f, x0, R0)
% maps - computes a Taylor series for mapping a set under linear dynamics
% onto a halfspace
%
% Syntax:  
%    maps(obj, A, f)
%
% Inputs:
%    obj - halfspace object
%    A - system matrix
%    f - contant input
%    x0 - center of initial states
%    x0int - interval vector of possible initial states
%
% Outputs:
%    k - constant map
%    L - linear map
%    Q - quadratic map
%    C - cubic map
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      14-June-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get normal vector and distance
n = obj.c;
d = obj.d;

%identity
I = eye(length(A));

%compute auxiliary variables
Lambda = n'*(A*x0 + f);
Upsilon = (n'*A)';
Theta = -n*Lambda - (d - n'*x0)*Upsilon;
Omega = -n*Upsilon' + Upsilon*n';

%compte constant vector k
k = x0 + (A*x0 + f)*(d - n'*x0)/Lambda;

%compte linear map L
L = I + ...
    A*(d-n'*x0)/Lambda + ...
    (A*x0 + f)*Theta'/Lambda^2;

%quadratic map Q(i,l,m)
for i=1:length(f)
    Q{i} =  (A(i,:)'*Theta' + Theta*A(i,:))/Lambda^2 ...
          + (A(i,:)*x0 + f(i))*(Omega*Lambda - Theta*2*Upsilon')/Lambda^3;
end

%auxiliary result
%aux = get(zonotope(interval(A*R0 + f)),'Z');
aux = get(A*R0 + f,'Z');

%cubic map C(i,l,m,n)
for i=1:length(f)
    for l=1:length(f)
        for m=1:length(f)
            for n=1:length(f)
                %compute center
                C_c{i}(l,m,n) =  (...
                               (A(i,l)*Omega(m,n) + A(i,m)*Omega(l,n))*Lambda ...
                             - (A(i,l)*Theta(m) + A(i,m)*Theta(l))*2*Upsilon(n) ...
                             + A(i,n)*(Omega(l,m)*Lambda - Theta(l)*2*Upsilon(m)) ...
                               )/Lambda^3 ...
                             + aux(i,1)*(...
                               (Omega(l,m)*Upsilon(n) - Omega(l,n)*2*Upsilon(m))*Lambda ...
                             - (Omega(l,m)*Lambda - Theta(l)*2*Upsilon(m))*3*Upsilon(n) ...
                               )/Lambda^4;
                 %compute delta  
                 for iGen = 1:length(aux(1,2:end))                     
                     C_g{i,iGen}(l,m,n) = aux(i,iGen+1)*(... 
                               (Omega(l,m)*Upsilon(n) - Omega(l,n)*2*Upsilon(m))*Lambda ...
                             - (Omega(l,m)*Lambda - Theta(l)*2*Upsilon(m))*3*Upsilon(n) ...
                               )/Lambda^4;                      
                 end
            end
        end
    end
end

%------------- END OF CODE --------------