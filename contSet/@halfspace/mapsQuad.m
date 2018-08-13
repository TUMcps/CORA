function [k,L,Q,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(obj, A, f, x0, R0, varargin)
% mapsQuad - computes a Taylor series for mapping a set under linear 
% dynamics onto a halfspace
%
% Syntax:  
%    [k,L,Q,lambdaQuad_inv_zono,Lambda_int] = mapsQuad(obj, A, f, x0, R0, varargin)
%
% Inputs:
%    obj - halfspace object
%    A - system matrix
%    f - contant input
%    x0 - center of initial states
%    R0 - set of initial states
%
% Outputs:
%    k - constant map
%    L - linear map
%    Q - quadratic map
%    lambdaQuad_inv_zono - zonotope of the inverse of the quadratic value
%    of Lambda
%    Lambda_int - interval value of Lambda
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
%               06-September-2013
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 5
    tilde_n = obj.c;
    tilde_d = obj.d;
else
    tilde_n = varargin{1};
    tilde_d = varargin{2};
end

%get normal vector of halfspace
n = obj.c;

%compute auxiliary variables
Lambda = n'*(A*x0 + f);
Upsilon = n'*A;
Theta = -tilde_n*Lambda - (tilde_d - tilde_n'*x0)*Upsilon';
Omega = -tilde_n*Upsilon + Upsilon'*tilde_n';

%set-based values
Z_mat = get(R0,'Z');
Lambda_zono = n'*(A*R0 + f);
Theta_aux_zono = (-1)*((-1)* tilde_n'*R0 + tilde_d)*Upsilon'; % order of summation important due to subsequent exact addition
Theta_zono = exactPlus(-tilde_n*Lambda_zono, Theta_aux_zono); % use with caution!
Theta_mat = get(Theta_zono,'Z');

%interval versions of auxiliary variables
Lambda_int = interval(Lambda_zono);

%check if Lambda_int contains 0; if yes, an alternative guard intersection
%method has to be called
if ~(infimum(Lambda_int)<0 && supremum(Lambda_int)>=0)
    Theta_int = interval(Theta_zono);

    %obtain lambdaQuad_zono
    lambdaQuad_mid = mid(1/Lambda_int^2);
    lambdaQuad_rad{1} = rad(1/Lambda_int^2);
    lambdaQuad_inv_zono = matZonotope(lambdaQuad_mid, lambdaQuad_rad);

    %interval of Theta/Lambda
    Theda_int = Theta_int/Lambda_int;
    Theda_mid = mid(Theda_int);
    Theda_rad = rad(Theda_int);

    %compte constant vector k
    k = (A*x0 + f)*(tilde_d - tilde_n'*x0)/Lambda; %changed: x0 outsourced

    %compte linear map L
    L = A*(tilde_d-tilde_n'*x0)/Lambda + ... %changed: I outsourced
        (A*x0 + f)*Theta'/Lambda^2;

    %quadratic map Q(i,l,m)
    c = center(A*R0 + f);
    Z_mat = get(A*R0 + f,'Z');
    G = Z_mat(:,2:end);
    gens = length(G(1,:));
    for i=1:length(f)
        %compute center matrix
        c_Q = A(i,:)'*Theta_mat(:,1)' + Theta_mat(:,1)*A(i,:) ...
              + c(i)*(Omega - Theda_mid*2*Upsilon);
        c_Q_rad = c(i)*(- Theda_rad*2*Upsilon);

        %compute generator matrices    
        for iGen = 1:gens
            g_Q{iGen} = A(i,:)'*Theta_mat(:,1+iGen)' + Theta_mat(:,1+iGen)*A(i,:) ...
              + G(i,iGen)*(Omega - Theda_mid*2*Upsilon);
          
            g_Q_rad{iGen} = G(i,iGen)*(- Theda_rad*2*Upsilon);
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
