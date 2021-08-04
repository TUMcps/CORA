function Zres = intersectStrip(obj,C,phi,y,varargin)
% intersectStrip - computes the intersection between a constrained zonotope 
% and a list of strips according to [1] ; a strip is defined as 
% | Cx-y | <= phi
%
% Syntax:  
%    Zres = intersectStrip(Z,C,phi,y,varargin)
%
% Inputs:
%    obj - constrained zonotope object
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
%    res - boolean whether obj is contained in Z, or not
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
%    poly = mptPolytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);
%    Zpoly = cZ & poly;
% 
%    figure; hold on 
%    plot(cZ,[1 2],'r-+');
%    plot(poly,[1 2],'r-*');
%    plot(Zpoly,[1 2],'b-+');
%    plot(res_zono,[1 2],'b-*');
% 
%    legend('conZonotope','strips','conZono&poly','zonoStrips');
%
%
% References:
%   [1] Amr Alanwar, Victor Gassmann, Xingkang He, Hazem Said,
%       Henrik Sandberg, Karl Henrik Johansson, and Matthias
%       Althoff. Privacy preserving set-based estimation using
%       partially homomorphic encryption. arXiV.org.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-Mar-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==4
    %The optimization function is based on norm of the generators
    method = 'normGen';
elseif nargin==5
    % input is Lambda value
    if isnumeric(varargin{1})
        Zres = conZonotopeFromLambda(obj,phi,C,y,varargin{1});
        return % stop program execution
    % input is a struct
    elseif isstruct(varargin{1})
        method = varargin{1}.method;
        aux = varargin{1};
    % input is intersection method
    else
        method = varargin{1};
    end
end

% extract generators
G = obj.Z(:,2:end);
% obtain dimesnion and nr of generators
[dim, gens] = size(G);

% different methods for finding good Lambda values
%% svd and radius method
if strcmp(method,'svd') || strcmp(method,'radius') 
    % initialize Lambda
    Lambda0 = zeros(dim,length(phi));
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find optimal Lambda
    Lambda = fminunc(@fun,Lambda0, options);
    % resulting constrained zonotope
    Zres = conZonotopeFromLambda(obj,phi,C,y,Lambda);
    
    
%% F-radius minimization according to Alamo, [3]    
elseif strcmp(method,'alamo-FRad') 
    
    % auxiliary variables
    aux1 = G*G'; 
    aux2 = aux1*C(1,:)';
    aux3 = C(1,:)*aux1*C(1,:)' + phi(1)^2; 

    % return Lambda
    Lambda = aux2/aux3;
    
    % resulting constrained zonotope
    Zres = conZonotopeFromLambda(obj,phi,C,y,Lambda);
    
    % warning
    if length(phi) > 1
        disp('Alamo method should only be used for single strips to ensure convergence');
    end
    
%% F-radius minimization according to Wang, Theorem 2 in [4]    
elseif strcmp(method,'wang-FRad') 
    
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
    Lambda = L*inv(S);
    
    % resulting constrained zonotope
    Zres = conZonotopeFromLambda(obj,phi,C,y,Lambda);
    
  
%% method according to Bravo, [5]   
elseif strcmp(method,'bravo') 
    
    % obtain center of constrained zonotope
    p = center(obj);
    % loop through generators
    for j = 0:gens
        % first generator 
        if j == 0
            v = p; % new center
            T = G; % new generators
        % normal vector of strip and generator are not perpendicular
        elseif  abs(C(1,:)*G(:,j)) > 10*eps
            % init T
            T = zeros(dim,gens);
            for i =1:gens
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
        volApprxZ(j+1)  = det(G*G');
    end
    % find zonotope with minimum approximated volume
    [~,ind] = min(volApprxZ);

    % return best zonotope as a constrained zonotope
    Zres = conZonotope(Z_candidate{ind});
    
    % warning
    if length(phi) > 1
        disp('Bravo method is only applied to the first strip');
    end


%% norm Gen method
elseif strcmp(method,'normGen')
    % Find the analytical solution  
    gamma=eye(length(C(:,1)));
    num= G*G'*C';
    den = C*G*G'*C';
    for iStrip=1:length(C(:,1))
        den = den + gamma(:,iStrip) *phi(iStrip)^2* gamma(:,iStrip)';
    end
    
    Lambda = num * den^-1;
    % resulting constrained zonotope
    Zres = conZonotopeFromLambda(obj,phi,C,y,Lambda);
    
    
%% selected method does not exist
else
    disp('Method is not supported');
    return;
end


    % embedded function to be minimized for optimal Lambda
    function nfro = fun(Lambda)
        part1 = eye(length(obj.center));
        for i=1:length(phi)
            part1 = part1 - Lambda(:,i)*C(i,:);
            part2(:,i) = phi(i)*Lambda(:,i);
        end
        part1 = part1 * G;
        G_new = [part1 part2];
        if strcmp(method,'svd')
            nfro = sum(svd(G_new));
        elseif strcmp(method,'radius')
            nfro = radius(zonotope([zeros(length(obj.center),1) G_new]));
        end

    end

end

% return constrained zonotope from a given Lambda vector, see Theorem 6.3. of [1]
function obj = conZonotopeFromLambda(obj,phi,C,y,Lambda)
    % strips: |Cx − y| <= phi
    
    % number of strips, dimensions
    nrOfStrips = length(phi);
    dim = length(obj.Z(:,1));
    
    % new center
    c = obj.Z(:,1);
    c_new = c + Lambda*(y-C*c);

    % compute new generators
    part1 = eye(dim);
    if isempty(obj.A)
        A_new =[];
        b_new =[];
    else
        A_new = [obj.A, zeros(size(obj.A,1),length(phi))];
        b_new = obj.b;
    end

    for i=1:nrOfStrips
        part1 = part1 - Lambda(:,i)*C(i,:);
        part2(:,i) = phi(i)*Lambda(:,i);
        A_new = [A_new ; C(i,:)*obj.Z(:,2:end) , zeros(1,i-1),-phi(i),zeros(1,length(phi)-i)];
        b_new = [b_new; y(i)-(C(i,:)*obj.Z(:,1))];
    end
    part1 = part1 * obj.Z(:,2:end);
    H_new = [part1 part2];

    % resulting constrained zonotope
    obj.Z = [c_new H_new];
    obj.A = A_new;
    obj.b = b_new;

end

%------------- END OF CODE --------------
