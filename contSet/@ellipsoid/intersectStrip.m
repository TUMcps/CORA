function [E,sigma_sq] = intersectStrip(E,C,phi,y,varargin)
% intersectStrip - computes the intersection between an ellipsoid and a 
%    list of strips, where a strip is defined as | Cx-y | <= phi
%
% Syntax:
%    Eres = intersectStrip(E,C,phi,y,varargin)
%
% Inputs:
%    E - ellipsoid object
%    C - matrix of normal vectors of strips
%    phi - vector of widths of strips/auxiliary struct for Liu2016
%    y - center of intersected strips
%    sigma_seq_prev - 
%    methods to calculate the weights
%               'Gollamudi1996' according to [1]
%               'Liu2016' according to [2]
%
% Outputs:
%    E - enclosing ellipsoid
%    sigma_sq
%
% Example: (one strip and one ellipsoid)
%    C = [1 0];
%    phi = 0.5;
%    y = 0;
%    sigma_sq_prev = 0.5;
% 
%    % convert strip to polytope
%    C_poly = [C; -C];
%    d = [phi + y; phi - y];
%    P = polytope(C_poly,d);
% 
%    E = ellipsoid([1 0.5; 0.5 1]);
%    Eres = intersectStrip(E,C,phi,y,sigma_sq_prev);
% 
%    figure; hold on 
%    plot(E,[1 2],'r');
%    plot(P,[1 2],'k');
%    plot(Eres,[1 2],'b');
%    legend('ellipsoid','strip','enclosing ellipsoid');
%
% References:
%    [1] S. Gollamudi, S. Nagaraj, S. Kapoor, and Y. F. Huang.
%        Set-membership state estimation with optimal bounding
%        ellipsoids. In Proc. of the International Symposium on
%        Information Theory and its Applications, pages 262–265,
%        1996.
%    [2] Yushuang Liu, Yan Zhao, and Falin Wu. Ellipsoidal state-
%        bounding-based set-membership estimation for linear system
%        with unknown-but-bounded disturbances.
%        IET Control Theory & Applications, 10(4):431–442, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-January-2021
% Last update:   04-July-2022 (VG, input checks)
%                29-March-2023 (TL, clean up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% parse input
if nargin < 4
    throw(CORAerror("CORA:notEnoughInputArgs", 4))
elseif nargin > 6
    throw(CORAerror('CORA:tooManyInputArgs', 6))
end
[sigma_sq_prev,method] = setDefaultValues({0,'Gollamudi1996'},varargin);
if ~isnumeric(sigma_sq_prev)
    method = sigma_sq_prev;
    sigma_sq_prev = 0;
end

inputArgsCheck({{E,'att','ellipsoid'}; ...
    {C, 'att', 'numeric', 'matrix'}; ...
    {phi, 'att', {'numeric','struct'}}; ...
    {y, 'att', 'numeric', 'column'}; ...
    {sigma_sq_prev, 'att' 'numeric', 'scalar'}; ...
    {method, 'str', {'Gollamudi1996', 'Liu2016'}}
})

% different methods for finding good lambda values
if strcmp(method,'Gollamudi1996')  
    if ~isnumeric(phi) && ~iscolumn(phi)
        throw(CORAerror("CORA:wrongValue", 'third', 'numeric column'))
    end

    [E,sigma_sq] = aux_methodGollamudi1996(E, C, phi, y, sigma_sq_prev);
    
elseif strcmp(method,'Liu2016')   
    if ~isstruct(phi)
        throw(CORAerror("CORA:wrongValue", 'third', 'struct'))
    end

    [E,sigma_sq] = aux_methodLiu2016(E, C, phi, y, sigma_sq_prev);
else
    throw(CORAerror('CORA:wrongValue', sprintf("Unknown method '%s'.", method)))
end

end


% Auxiliary functions -----------------------------------------------------

function [E,sigma_sq] = aux_methodGollamudi1996(E, C, phi, y, sigma_sq_prev)
    % intersection according to Gollamudi1996, [1]
    % ||Cx-y|| <= phi
    % correspondences to [1]:
    % \gamma --> phi
    % x --> c (center)
    % y is identical
    % P --> Q (shape matrix)

    % obtain center and shape matrix
    c = center(E);
    Q = E.Q;
    
    % auxiliary value
    delta = y - C*c;
    
    % identity matrix
    I = eye(length(y));

    % Theorem 2 of [1]
    if sigma_sq_prev + delta'*delta > phi

        % fix alpha value
        alpha = 1e-4;

        % G, see Theorem 1 of [1]
        G = C*Q*C';

        % maximum singular value of G
        g = svds(G,1);

        % compute nu of Theorem 2 in [1]
        % order is important: first check delta == 0
        if withinTol(delta,0)
            nu = alpha;
        else
            % beta, eq. (18) of [1]
            beta = (phi - sigma_sq_prev)/(delta'*delta);
            if withinTol(g,1)
                nu = (1-beta)/2;
            else
                %compute omega
                omega = 1 + beta*(g-1);
                if omega > 0
                    nu = 1/(1-g)*(1-sqrt(g/omega));
                else
                    nu = alpha;
                end
            end
        end

        % obtain lambda
        lambda = min(alpha, nu);

        % new center
        % original: x := x + lambda*P*C'*delta
        c = c + lambda*Q*C'*delta;

        % new \sigma^2 used for scaling; sigma=1 here
        % original: \sigma^2 = (1-\lambda)*\sigma^2 + \lambda \gamma^2 -
        % \lambda(1-\lambda)\delta'*G^(-1)*delta
        % original: G = (1-\lambda)*I + \lambda*C*P*C'
        Q_tmp = (1-lambda)*I + lambda*G;
        sigma_sq = 1-lambda + lambda*phi^2 - lambda*(1-lambda)*delta'*inv(Q_tmp)*delta;

        % new shape matrix
        % original: inv(P) = (1-\lambda)inv(P) + \lambda*C'*C
        invQ = (1-lambda)*inv(Q) + lambda*(C'*C);
        Q = inv(invQ)/sigma_sq;

        % resulting ellipsoid
        E.Q = Q;
        E.q = c;
    else
        % ellipsoid is unchanged
        sigma_sq = sigma_sq_prev;
    end
end

function [E,sigma_sq] = aux_methodLiu2016(E, C, sys, y, sigma_sq_prev)
    % intersection according to Liu2016, [2]    
    % correspondences to [2]:
    % \gamma --> phi
    % x --> c (center)
    % y is identical
    % P --> Q (shape matrix)
    % F --> A
    % H --> C

    % obtain center and shape matrix
    c = center(E);
    Q = E.Q;
    
    % auxiliary value
    delta = y - C*c;

    % Theorem 3 of [2]
    if sigma_sq_prev + delta'*delta > 1

        % G, see below eq. (51) of [2]
        G = sys.bar_V*C*Q*C'*sys.bar_V';

        % maximum eigenvalue of G
        g = eigs(G,1);

        % compute lambda of Theorem 3 in [2]
        % beta, eq. (52) of [2]
        beta = (1 - sigma_sq_prev)/(delta'*delta);
        if withinTol(g,1)
            %compute lambda
            lambda = (1-beta)/2;
        else
            %compute lambda
            omega = 1 + beta*(g-1);
            lambda = 1/(1-g)*(1-sqrt(g/omega));
        end
        
        % eq. (27) of [2]
        Q_tmp = 1/lambda*sys.V + 1/(1-lambda)*C*Q*C';
        
        % eq. (26) of [2]
        Q_tmp_inv = inv(Q_tmp);
        K = 1/(1-lambda)*Q*C'*Q_tmp_inv;
        
        % eq. (25) of [2]
        sigma_sq = (1-lambda)*sigma_sq_prev + lambda - delta'*Q_tmp_inv*delta;
            
        % new shape matrix, eq. (24) of [2]
        I = eye(length(Q)); % identity matrix
        Q = (1-lambda)*(I - K*C)*Q;

        % new center, eq. (23) of [2]
        c = c + K*delta;

        % resulting ellipsoid
        E.Q = Q;
        E.q = c;
    else
        % ellipsoid is unchanged
        sigma_sq = sigma_sq_prev;
    end
end

% ------------------------------ END OF CODE ------------------------------
