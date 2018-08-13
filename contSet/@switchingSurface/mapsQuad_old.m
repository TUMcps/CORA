function [k,L,Q,C] = mapsQuad_old(obj, A, f, x0, R0)
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
aux = get(A*R0 + f,'Z');
for i=1:length(f)
    %compute center matrix
    c_Q = (A(i,:)'*Theta' + Theta*A(i,:))/Lambda^2 ...
          + aux(i,1)*(Omega*Lambda - Theta*2*Upsilon')/Lambda^3;
      
    %compute generator matrices    
    for iGen = 1:length(aux(1,2:end))
        g_Q{iGen} =  aux(i,1+iGen)*(Omega*Lambda - Theta*2*Upsilon')/Lambda^3;
    end
      
    %generate matrix zonotope
    Q{i} = matZonotope(c_Q, g_Q);
end


%------------- END OF CODE --------------