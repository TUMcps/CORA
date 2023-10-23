function cZ = intersectStrip(cZ,C,phi,y,varargin)
% intersectStrip - computes the intersection between a constrained zonotope 
%    and a list of strips according to [1] ; a strip is defined as 
%    | Cx-y | <= phi
%
% Syntax:
%    Zres = intersectStrip(Z,C,phi,y,varargin)
%
% Inputs:
%    cZ - conZonotope object
%    C - matrix of normal vectors of strips
%    phi - vector of widths of strips
%    y - center of intersected strips
%    varargin - methods to calculate the weights
%               'normGen' default and has analytical solution
%               'svd'
%               'radius'
%               'alamo-FRad' according to [3]
%               'wang-FRad' according to [4] + auxiliary values as a struct
%               'bravo' accroding to [5]
%               or Lambda value
%
% Outputs:
%    cZ - conZonotope object
%
% Example: (three strips and one constrained zonotope)
%    C = [1 0; 0 1; 1 1];
%    phi = [5; 3; 3];
%    y = [-2; 2; 2];
% 
%    Z = [1 2 2 2 6 2 8; 1 2 2 0 5 0 6];
%    A = [1 0 1 0 0 0];
%    b = 1;
%    cZ = conZonotope(Z,A,b);
%    res_zono = intersectStrip(cZ,C,phi,y);
% 
%    % just for comparison
%    poly = polytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);
%    Zpoly = cZ & poly;
% 
%    figure; hold on;
%    plot(cZ,[1 2],'b','DisplayName','conZonotope');
%    plot(poly,[1 2],'k*','DisplayName','Strips');
%    plot(Zpoly,[1 2],'r','DisplayName','conZono&strip');
%    plot(res_zono,[1 2],'g--','DisplayName','zonoStrips');
%    legend()
%
%
% References:
%   [1] Amr Alanwar, Victor Gassmann, Xingkang He, Hazem Said,
%       Henrik Sandberg, Karl Henrik Johansson, and Matthias
%       Althoff. Privacy preserving set-based estimation using
%       partially homomorphic encryption. arXiv.org.
%   [2-5] ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-March-2020
% Last update:   29-March-2023 (TL, clean up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 4
    throw(CORAerror("CORA:notEnoughInputArgs", 4))
elseif nargin > 5
    throw(CORAerror("CORA:tooManyInputArgs", 5))
end

inputArgs = { ...
    {cZ, 'att', 'conZonotope'}; ...
    {C, 'att', 'numeric', 'matrix'}; ...
    {phi, 'att', 'numeric', 'column'}; ...
    {y, 'att', 'numeric', 'column'}; ...
    };
inputArgsCheck(inputArgs)

% parse varargin 
method = setDefaultValues({'normGen'}, varargin);
if isnumeric(method)
    lambda = method;
    cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,lambda);
    return % stop program execution

elseif isstruct(method)
    aux = method;
    aux_checkAuxStruct(aux)
    method = aux.method;
end

inputArgsCheck([inputArgs; ...
    {{method, 'str', {'normGen', 'svd', 'radius', 'alamo-FRad', ...
    'wang-FRad', 'bravo'}}} ...
])

% check dimensions
if CHECKS_ENABLED
    if dim(cZ) ~= size(C, 2)
        throw(CORAerror('CORA:dimensionMismatch', cZ, C))
    end
    if size(C, 1) ~= size(phi, 1)
        throw(CORAerror('CORA:dimensionMismatch', C, phi))
    end
    if size(C, 1) ~= size(y, 1)
        throw(CORAerror('CORA:dimensionMismatch', C, y))
    end
end

% different methods for finding good Lambda values
if strcmp(method,'svd') || strcmp(method,'radius')
    cZ = aux_methodSvd(cZ, C, phi, y, method);

elseif strcmp(method,'alamo-FRad') 
    cZ = aux_methodAlamoFRad(cZ, C, phi, y);
    
elseif strcmp(method,'wang-FRad') 
    cZ = aux_methodWangFRad(cZ, C, phi, y, aux);
    
elseif strcmp(method,'bravo') 
    cZ = aux_methodBravo(cZ, C, phi, y);

elseif strcmp(method,'normGen')
    cZ = aux_methodNormGen(cZ, C, phi, y);

else
    throw(CORAerror('CORA:wrongValue', sprintf("Unknown method '%s'.", method)))
end

end


% Auxiliary functions -----------------------------------------------------

% method functions ---

function cZ = aux_methodSvd(cZ, C, phi, y, method)
    % svd and radius method

    G = generators(cZ);
    [n, ~] = size(G);

    % initialize Lambda
    Lambda0 = zeros(n,length(phi));
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find optimal Lambda
    Lambda = fminunc(@fun,Lambda0, options);
    % resulting constrained zonotope
    cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,Lambda);

    function nfro = fun(Lambda)
        % embedded function to be minimized for optimal Lambda
        % (used for svd/radius method)
        part1 = eye(length(cZ.center));
        for i=1:length(phi)
            part1 = part1 - Lambda(:,i)*C(i,:);
            part2(:,i) = phi(i)*Lambda(:,i);
        end
        part1 = part1 * G;
        G_new = [part1 part2];
        if strcmp(method,'svd')
            nfro = sum(svd(G_new));
        elseif strcmp(method,'radius')
            nfro = radius(zonotope([zeros(length(cZ.center),1) G_new]));
        end
    end
end

function cZ = aux_methodAlamoFRad(cZ, C, phi, y)
    % F-radius minimization according to Alamo, [3]

    G = generators(cZ);
    
    % auxiliary variables
    aux1 = G*G'; 
    aux2 = aux1*C(1,:)';
    aux3 = C(1,:)*aux1*C(1,:)' + phi(1)^2; 

    % return Lambda
    Lambda = aux2/aux3;
    
    % resulting constrained zonotope
    cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,Lambda);
    
    % warning
    if length(phi) > 1
        warning('Alamo method should only be used for single strips to ensure convergence');
    end
end

function cZ = aux_methodWangFRad(cZ, C, phi, y, aux)
    % F-radius minimization according to Wang, Theorem 2 in [4] 

    G = generators(cZ);
    
    % auxiliary variables
    P = G*G'; 
    Q_w = aux.E*aux.E';
    Q_v = aux.F*aux.F'; 

    % eq. (15)
    Rbar = aux.A*P*aux.A' + Q_w;
    
    % eq. (14)
    S = aux.C*Rbar*aux.C' + Q_v;
    
    % eq. (13)
    L = Rbar*aux.C';
    
    % eq. (12)
    Lambda = L / S; %L*inv(S);
    
    % resulting constrained zonotope
    cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,Lambda);
end


function cZ = aux_methodBravo(cZ, C, phi, y)
    % method according to Bravo, [5]   

    G = generators(cZ);
    % obtain dimesnion and nr of generators
    [n, nrGens] = size(G);

    % obtain center of constrained zonotope
    p = center(cZ);

    Z_candidate = cell(1, nrGens+1);
    volApprxZ = zeros(1, nrGens+1);

    % loop through generators
    for j = 0:nrGens
        % first generator 
        if j == 0
            v = p; % new center
            T = G; % new generators
        % normal vector of strip and generator are not perpendicular
        elseif  abs(C(1,:)*G(:,j)) > 10*eps
            % init T
            T = zeros(n,nrGens);
            for i =1:nrGens
                if i~=j
                    T(:,i) = G(:,i)  - C(1,:)*G(:,i)/(C(1,:)*G(:,j))*G(:,j);
                elseif i == j
                    T(:,i) = phi(1)/(C(1,:)*G(:,j)).*G(:,j);
                end
            end
            v =  p + (y(1)-C(1,:)*p)./(C(1,:)*G(:,j))*G(:,j);
        % normal vector of strip and generator are perpendicular
        else
            v = p; % new center
            T = G; % new generators
        end
        % save new zonotope
        Z_candidate{j+1} = zonotope([v,T]);
        % approximate volume of obtained zonotope
        G = generators(Z_candidate{j+1});
        volApprxZ(j+1) = det(G*G');
    end
    % find zonotope with minimum approximated volume
    [~,ind] = min(volApprxZ);

    % return best zonotope as a constrained zonotope
    cZ = conZonotope(Z_candidate{ind});
    
    % warning
    if length(phi) > 1
        warning('Bravo method is only applied to the first strip');
    end
end

function cZ = aux_methodNormGen(cZ, C, phi, y)
    % norm Gen method
        
     G = generators(cZ);

    % Find the analytical solution  
    gamma=eye(length(C(:,1)));
    num= G*G'*C';
    den = C*G*G'*C';
    for iStrip=1:length(C(:,1))
        den = den + gamma(:,iStrip) *phi(iStrip)^2* gamma(:,iStrip)';
    end
    
    Lambda = num * den^-1;
    % resulting constrained zonotope
    cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,Lambda);
end

% helper functions ---

% return constrained zonotope from a given Lambda vector, see Theorem 6.3. of [1]
function cZ = aux_conZonotopeFromLambda(cZ,phi,C,y,Lambda)
    % strips: |Cx − y| <= phi
    
    % number of strips, dimensions
    nrOfStrips = length(phi);
    n = dim(cZ);
    
    % new center
    c = cZ.c;
    c_new = c + Lambda*(y-C*c);

    % compute new generators
    part1 = eye(n);
    if isempty(cZ.A)
        A_new =[];
        b_new =[];
    else
        A_new = [cZ.A, zeros(size(cZ.A,1),length(phi))];
        b_new = cZ.b;
    end

    for i=1:nrOfStrips
        part1 = part1 - Lambda(:,i)*C(i,:);
        part2(:,i) = phi(i)*Lambda(:,i);
        A_new = [A_new; C(i,:)*cZ.G , zeros(1,i-1),-phi(i),zeros(1,length(phi)-i)];
        b_new = [b_new; y(i)-(C(i,:)*cZ.c)];
    end
    part1 = part1 * cZ.G;
    H_new = [part1 part2];

    % resulting constrained zonotope
    cZ.c = c_new;
    cZ.G = H_new;
    cZ.A = A_new;
    cZ.b = b_new;

end

function aux_checkAuxStruct(aux)
    if CHECKS_ENABLED
        if ~isfield(aux, 'method')
            throw(CORAerror('CORA:wrongFieldValue', 'aux.method', {'wang-FRad'}))
        end
        if ~isfield(aux, 'E')
            throw(CORAerror('CORA:wrongFieldValue', 'aux.E', 'numeric'))
        end
        if ~isfield(aux, 'F')
            throw(CORAerror('CORA:wrongFieldValue', 'aux.F', 'numeric'))
        end
        if ~isfield(aux, 'A')
            throw(CORAerror('CORA:wrongFieldValue', 'aux.A', 'numeric'))
        end
        if ~isfield(aux, 'C')
            throw(CORAerror('CORA:wrongFieldValue', 'aux.C', 'numeric'))
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
