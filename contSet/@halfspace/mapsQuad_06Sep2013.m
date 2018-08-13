function [k,L,Q,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(obj, A, f, x0, R0)
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
% Last update:  28-August-2013
%               25-July-2016 (intervalhull replaced by interval)
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

%set-based values
Lambda_zono = n'*(A*R0 + f);
Theta_aux_zono = (-1)*(d + (-1)* n'*R0)*Upsilon;
Theta_aux_mat = get(Theta_aux_zono,'Z');
Lambda_mat = get(Lambda_zono,'Z');
Theta_mat = -n*Lambda_mat+Theta_aux_mat;

%interval versions of auxiliary variables
Lambda_int = interval(Lambda_zono);

%check if Lambda_int contains 0; if yes, an alternative guard intersection
%method has to be called
if ~(infimum(Lambda_int)<0 && supremum(Lambda_int)>=0)
    Theta_int = interval(zonotope(Theta_mat));

    %obtain lambdaQuad_zono
    lambdaQuad_mid = mid(1/Lambda_int^2);
    lambdaQuad_rad{1} = rad(1/Lambda_int^2);
    lambdaQuad_inv_zono = matZonotope(lambdaQuad_mid, lambdaQuad_rad);

    %interval of Theta/Lambda
    Theda_int = Theta_int/Lambda_int;
    Theda_mid = mid(Theda_int);
    Theda_rad = rad(Theda_int);

    %compte constant vector k
    k = x0 + (A*x0 + f)*(d - n'*x0)/Lambda;

    %compte linear map L
    L = I + ...
        A*(d-n'*x0)/Lambda + ...
        (A*x0 + f)*Theta'/Lambda^2;

    %quadratic map Q(i,l,m)
    c = center(A*R0 + f);
    Z_mat = get(A*R0 + f,'Z');
    G = Z_mat(:,2:end);
    gens = length(G(1,:));
    for i=1:length(f)
        %compute center matrix
        c_Q = A(i,:)'*Theta_mat(:,1)' + Theta_mat(:,1)*A(i,:) ...
              + c(i)*(Omega - Theda_mid*2*Upsilon');
        c_Q_rad = c(i)*(- Theda_rad*2*Upsilon');

        %compute generator matrices    
        for iGen = 1:gens
            g_Q{iGen} = A(i,:)'*Theta_mat(:,1+iGen)' + Theta_mat(:,1+iGen)*A(i,:) ...
              + G(i,iGen)*(Omega - Theda_mid*2*Upsilon');
          
            g_Q_rad{iGen} = G(i,iGen)*(- Theda_rad*2*Upsilon');
        end

        Q_prep{1} = c_Q_rad;
        Q_prep(2:(gens+1)) = g_Q;
        Q_prep((gens+2):(2*gens+1)) = g_Q_rad;

        %generate matrix zonotope
        Q{i} = matZonotope(c_Q, Q_prep);
    end
else
    k = [];
    L = [];
    Q = [];
    lambdaQuad_inv_zono = [];
    Lambda_int = [];
end

%------------- END OF CODE --------------