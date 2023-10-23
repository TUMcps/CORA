function cZ = quadMap(cZ,varargin)
% quadMap - computes the quadratic map of a constrained zonotope
%
% Syntax:
%    cZquad = quadMap(cZ,Q)
%    cZquad = quadMap(cZ,cZ2,Q)
%
% Inputs:
%    cZ - conZonotope object
%    cZ2 - conZonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    Q{1} = [1 2;-1 0];
%    Q{2} = [-2 -1;0 1];
%    cZquad = quadMap(cZ,Q);
%
%    figure; hold on;
%    plot(cZ,[1,2],'FaceColor','r');
%    plot(cZquad,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap

% Authors:       Niklas Kochdumper
% Written:       13-August-2018
% Last update:   22-November-2019 (NK, integrated mixed quad. mul.)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if nargin == 1
        throw(CORAerror('CORAerror:notEnoughInputArgs',2));
    elseif nargin == 2
        cZ = aux_quadMapSingle(cZ,varargin{1});
    elseif nargin == 3
        cZ = aux_quadMapMixed(cZ,varargin{1},varargin{2});
    else
        throw(CORAerror('CORAerror:tooManyInputArgs',3));
    end

end


% Auxiliary functions -----------------------------------------------------

function cZquad = aux_quadMapSingle(cZ,Q)
% comptue the quadratic map {x_i = x^T Q x | x \in Z} of a zonotope

    % rescale the constrained zonotope to reduce the over-approximation of
    % the quadratic map
    if ~isempty(cZ.A)
        cZ = rescale(cZ,'iter');
    end

    % get generator matrix of constrained zonotope
    Zmat = [cZ.c,cZ.G];
    dimQ = length(Q);
    gens = length(Zmat(1,:)) - 1;

    % initialize variables
    N = (gens+1) * 0.5*(gens+2)-1;
    c = zeros(dimQ,1);
    G = zeros(dimQ,N);

    % loop over all dimensions
    for i = 1:dimQ

        % quadratic evaluation
        qMat = Zmat'*Q{i}*Zmat;

        % center
        c(i,1) = qMat(1,1) + 0.5*sum(diag(qMat(2:end,2:end)));

        % generators from diagonal elements
        ind = 1:gens;
        G(i, ind) = 0.5*diag(qMat(ind+1,ind+1));

        % generators from other elements
        counter = 0;
        for j = 0:gens
            kInd = (j+1):gens;
            G(i, gens + counter + kInd - j) = qMat(j+1, kInd+1) + ...
                                              qMat(kInd+1, j+1)';
            counter = counter + length(kInd);
        end
    end

    % construct the constrained matrices for the constrained zonotope that
    % represents the resulting set for the quadratic map
    N_ = N - 2*gens;
    A = [zeros(size(cZ.A,1),gens), cZ.A, zeros(size(cZ.A,1),N_)];

    % generate the resulting conZonotope object
    cZquad = conZonotope([c, G],A,cZ.b);
end

function cZquad = aux_quadMapMixed(cZ1,cZ2,Q)
% comptue the quadratic map {x_i = x1^T Q x2 | x1 \in Z1, x2 \in Z2} of two
% zonotope objects

    % rescale the constrained zonotopes to reduce the over-approximation of
    % the quadratic map
    if ~isempty(cZ1.A)
        cZ1 = rescale(cZ1,'iter');
    end

    if ~isempty(cZ2.A)
        cZ2 = rescale(cZ2,'iter');
    end

    % get generator matrices of the constrained zonotopes
    Zmat1 = [cZ1.c,cZ1.G];
    Zmat2 = [cZ2.c,cZ2.G];
    dimQ = length(Q);

    % initialize variables
    N1 = size(Zmat1,2)-1;
    N2 = size(Zmat2,2)-1;
    c = zeros(dimQ,1);
    G = zeros(dimQ,(N1+1)*(N2+1)-1);

    % loop over all dimensions
    for i = 1:dimQ

        % quadratic evaluation
        quadMat = Zmat1'*Q{i}*Zmat2;

        % center 
        c(i,1) = quadMat(1,1);

        % generators resulting from a multiplication with a zonotope center
        G(i,1:N2) = quadMat(1,2:end);
        G(i,N2+1:N1+N2) = quadMat(2:end,1);

        % transform quadratic matrix to vector
        quadVec = reshape(quadMat(2:end,2:end),1,[]);

        % remaining generators
        G(i, N2+N1+1:end) = quadVec;
    end

    % construct the new constraint matrix   
    if isempty(cZ1.A) && isempty(cZ2.A)
        A = []; b = [];
    else
        A = [cZ2.A, zeros(size(cZ2.A,1),(N1+1)*(N2+1)-N2-1);
             zeros(size(cZ1.A,1),N2), cZ1.A, ...
             zeros(size(cZ1.A,1),(N1+1)*(N2+1)-N2-N1-1)];

        b = [cZ2.b;cZ1.b];
    end

    % generate the resulting constrained zonotope
    cZquad = conZonotope(c,G,A,b);

end

% ------------------------------ END OF CODE ------------------------------
