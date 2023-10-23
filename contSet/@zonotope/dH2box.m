function dH = dH2box(Z,varargin)
% dH2box - computes an over-approximation of the Hausdorff distance
%    to the interval over-approximation of the provided zonotope Z
%    note: this function is not optimized w.r.t computational efficiency
%
% Syntax:
%    dH = dH2box(Z,method)
%
% Inputs:
%    Z - zonotope object
%    method - char array of over-approximation method:
%       - 'exact' [eq.(5), 1]
%       - 'naive': radius of box over-approximation
%       - 'ell': sequence of parallelotopes P_i, each with inscribed
%                   ball containing box (= under-approximation of P_i)
%       - 'wgreedy' [Thm.3.2, 2]
%       - 'wopt': modification of [Thm.3.2, 2]
%
% Outputs:
%    dH - over-approximated Hausdorff distance
%
% Example:
%    Z = zonotope([0;0],-1+2*rand(2,20));
%    dH = dH2box(Z,'naive');
%
% References:
%    [1] X. Yang, J.K. Scott. A comparison of zonotope order reduction
%        techniques. Automatica 95 (2018), pp. 378-384.
%    [2] M. Wetzlinger, A. Kulmburg, M. Althoff. Adaptive parameter tuning
%        for reachability analysis of nonlinear systems. HSCC 2021.
%    [3] V. Gassmann, M. Althoff. Scalable zonotope-ellipsoid conversions
%        using the euclidean zonotope norm. ACC 2020.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Mark Wetzlinger
% Written:       08-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input arguments
if dim(Z) <= 3
    dV = 'exact';
else
    dV = 'naive';
end

% parse input arguments
method = setDefaultValues({dV},varargin);

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {method,'str',{'exact','naive','ell','wgreedy','wopt'}}});

% methods
if strcmp(method,'exact')
    % ... computes near-exact dH according to [1]
    dH = aux_dH2box_exact(Z);
elseif strcmp(method,'naive')
    % ... computes the radius of the box over-approximation of Z
    dH = vecnorm(rad(interval(Z)));
elseif strcmp(method,'ell')
    % ... computes the distance from a sequence of box under-approximations
    %       to the respective box over-approximations, based on ellipsoid
    %       inscriptions in zonotopes from [3]
    dH = aux_dH2box_ell(Z);
elseif strcmp(method,'wgreedy')
    % ... computes the length of the sum of the generators,
    %       where the largest entry (by absolute value) is set to 0 [2]
    dH = aux_dH2box_wgreedy(Z);
elseif strcmp(method,'wopt')
    % ... computes the length of the sum of the generators,
    %       where the length are scaled by choosing optimal factors
    dH = aux_dH2box_wopt(Z);
end


end


% Auxiliary functions -----------------------------------------------------

function dH = aux_dH2box_exact(Z)

    % generator matrices
    G = generators(Z);
    Gbox = diag(sum(abs(G),2));

    % generate lambda vectors
    n = dim(Z);
    numDirs = 10000;
    sections = nthroot(numDirs,n-1);
    lambda = randEqdistDirections(n,sections);
    
    % loop over each lambda to find maximum
    dH = 0;
    for i=1:size(lambda,2)
        dHlambda = abs(vecnorm(lambda(:,i)'*Gbox,1) - vecnorm(lambda(:,i)'*G,1));
        if dHlambda > dH
            dH = dHlambda;
        end
    end
    
    % caution! since this is an iterative max procedure, the result
    % is actually a (tight) under-approximation

end

function dH = aux_dH2box_ell(Z)

% generator matrix
G = generators(Z);
[n,gamma] = size(G);
% Gabs = abs(G);

% number of parallelotopes
nrP = floor(gamma/n);

% sort generators such that parallelotopes similar to boxes
% [maxval,maxdim] = max(Gabs,[],1);
% sumval = sum(Gabs,1);
% hgirard = sumval - maxval;

% take nrP*n longest generators
% Glength = vecnorm(G,2);
% [~,longestGens] = maxk(Glength,nrP*n);
% GforP = G(:,longestGens);
GforP = G;

% container for Hausdorff distances
dH = zeros(nrP,1);

% note: the function below is not optimized w.r.t computational efficiency
for i=1:nrP
    
    % extract invertible parallelotope P
    [~,p] = rref(GforP);
    if length(p) ~= n
        break;
    end
    P = GforP(:,p);
    % remove generators from G for parallelotopes
    GforP(:,p) = [];
    
    % vertex of interval over-approximation of P
    vPover = sum(abs(P),2);
    
    % invert generator matrix: A*x <= 1 is h-rep of P
    A = inv(P);
    
    % normalize rows of A and b
    c = 1./sqrt(sum(A.^2,2));
%     An = c.*A; % normalized rows of A not required
    bnorm = c;
    
    % radius of inscribed ball
    rball = min(abs(bnorm));
    
    % radius of box in inscribed ball
    rBox = sqrt(rball^2/n);
    
    % vertex of interval under-approximation of P
    vPunder = rBox*ones(n,1);
    
    % measure Hausdorff distance
    dH(i) = norm(vPover-vPunder,2);
    
    plotting = false;
    
    if plotting
        % plotting for visual understanding
        Z = zonotope(zeros(n,1),G);         % full zonotope
        Pt = zonotope(zeros(n,1),P);        % extracted parallelotope
        Pover = interval(Pt);               % over-approx of Pt
        E = ellipsoid(rball^2*eye(n));      % ellipsoid in Pt
        Punder = interval(-rBox*ones(n,1),rBox*ones(n,1)); % under-approx of Pt

        projDims = [1,2;3,4];
        rows = 1; cols = 2;
        figure;
        for proj = 1:size(projDims,1)
            subplot(rows,cols,proj); hold on; box on;
            currProj = projDims(proj,:);
            plot(Z,currProj,'k');
            plot(Pt,currProj,'b');
            plot(Pover,currProj,'g');
            plot(E,currProj,'r');
            plot(Punder,currProj,'r');
            xlabel("x_" + currProj(1));
            ylabel("x_" + currProj(2));
        end
        close;
    end
end

% accumulated Hausdorff distance
dH = sum(dH);

end

function dH = aux_dH2box_wgreedy(Z)

    % generator matrices
    G = generators(Z);
    Gabs = abs(G);
    [~,gamma] = size(G);

    % truncated generators (largest absolute element per column = 0)
    Gtruncated = Gabs;
    for k=1:gamma
        [~,idxmax] = max(Gabs(:,k));
        Gtruncated(idxmax,k) = 0;
    end
    
    % measure
    dH = 2*vecnorm(sum(Gtruncated,2),2);

end

function dH = aux_dH2box_wopt(Z)
    
    % generator matrices
    G = generators(Z);
    Gabs = abs(G);
    
    % dimensions
    [n,gamma] = size(G);
    
    % lengths of generators
    lsquared = vecnorm(G,2) .^ 2;
    
    dH = 0;
    for i=1:n
        summand = 0;
        for k=1:gamma
            if lsquared(k) ~= 0
                summand = summand + Gabs(i,k) * (1 - G(i,k)^2 / lsquared(k));
            end
        end
        dH = dH + summand^2;
    end
    dH = 2*sqrt(dH);
    
end

% ------------------------------ END OF CODE ------------------------------
