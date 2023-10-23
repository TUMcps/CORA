function Zquad = quadMap(Z1,varargin)
% quadMap - computes the quadratic map of a zonotope
%
% Syntax:
%    Zquad = quadMap(Z1,Q)
%    Zquad = quadMap(Z1,Z2,Q)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    Zquad - zonotope object
%
% Example: 
%    Z = zonotope([0 1 1;0 1 0]);
%    Q{1} = [0.5 0.5; 0 -0.5];
%    Q{2} = [-1 0; 1 1];
% 
%    res = quadMap(Z,Q);
% 
%    figure; hold on;
%    plot(Z,[1,2],'r');
% 
%    figure; hold on;
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems", HSCC 2011.
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       07-December-2011
% Last update:   22-November-2019 (NK, combined with mixed quad. mul.)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if nargin == 1
        throw(CORAerror('CORAerror:notEnoughInputArgs',2));
    elseif nargin == 2
        if isa(varargin{1}{1},'matZonotope')
            Zquad = aux_quadMapSingleMatZono(Z1,varargin{1});
        else
            Zquad = aux_quadMapSingle(Z1,varargin{1});
        end
    elseif nargin == 3
        Zquad = aux_quadMapMixed(Z1,varargin{1},varargin{2});
    else
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end
end


% Auxiliary functions -----------------------------------------------------

function Zquad = aux_quadMapSingle(Z,Q)
% compute an over-approximation of the quadratic map 
%
% {x_i = x^T Q{i} x | x \in Z} 
%
% of a zonotope according to Lemma 1 in [1]

    % get matrix of zonotope
    Zmat = [Z.c,Z.G];
    dimQ = length(Q);
    gens = size(Z.G,2);

    % init solution
    c = zeros(dimQ,1);
    G = zeros(dimQ,0.5*(gens^2+gens) + gens);
    
    % count empty matrices
    Qnonempty = false(dimQ,1);
    
    % for each dimension, compute generator elements
    for i = 1:dimQ

        Qnonempty(i) = any(Q{i}(:));
        if Qnonempty(i)
            
            % pure quadratic evaluation
            quadMat = Zmat'*Q{i}*Zmat;

            % faster method
            % diag elements
            G(i,1:gens) = 0.5*diag(quadMat(2:gens+1,2:gens+1));
            
            % center
            c(i,1) = quadMat(1,1) + sum(G(i,1:gens));

            % off-diagonal elements added, pick via logical indexing
            quadMatoffdiag = quadMat + quadMat';
            quadMatoffdiag = quadMatoffdiag(:);
            kInd = tril(true(gens+1,gens+1),-1);
            G(i, gens+1:end) = quadMatoffdiag(kInd(:));
        end

    end

    % generate new zonotope
    tmp_sum = sum(Qnonempty);
    if tmp_sum < 1 || withinTol(tmp_sum,1)  %if sum(Qnonempty) <= 1
        Zquad = zonotope([c,sum(abs(G),2)]);
    else
        Zquad = zonotope([c,nonzeroFilter(G)]);
    end
    
end

function Zquad = aux_quadMapMixed(Z1,Z2,Q)
% compute an over-approximation of the quadratic map 
%
% {x_i = x1^T Q{i} x2 | x1 \in Z1, x2 \in Z2} 
%
% of two zonotope objects.

    % get matrix of zonotope
    Zmat1 = [Z1.c,Z1.G];
    Zmat2 = [Z2.c,Z2.G];
    dimQ = length(Q);
    
    % init solution (center + generator matrix)
    Z = zeros(dimQ,size(Zmat1,2)*size(Zmat2,2));
    
    % count empty matrices
    Qnonempty = false(dimQ,1);

    % for each dimension, compute center + generator elements
    for i = 1:dimQ

        Qnonempty(i) = any(Q{i}(:));
        if Qnonempty(i)
            % pure quadratic evaluation
            quadMat = Zmat1'*Q{i}*Zmat2;
            Z(i,:) = quadMat(:)';
        end
        
    end

    % generate new zonotope
    tmp_sum = sum(Qnonempty);
    if tmp_sum < 1 || withinTol(tmp_sum,1)  %if sum(Qnonempty) <= 1
        Zquad = zonotope([Z(:,1),sum(abs(Z(:,2:end)),2)]);
    else
        Zquad = zonotope([Z(:,1),nonzeroFilter(Z(:,2:end))]);
    end

end

function Zquad = aux_quadMapSingleMatZono(Z,Q)
% compute an over-approximation of the quadratic map
%
% {x_i = x^T Q{i} x | x \in Z} 
%
% of a zonotope according to Theorem 1 in [1], where Q is a matrix zonotope

    % zonotope Z_D for the center of the matrix zonotope
    Q_ = cell(length(Q),1);
    
    for i = 1:length(Q)
       Q_{i} = Q{i}.center; 
    end
    
    Z_D = quadMap(Z,Q_);
    
    % zonotopes Z_Kj for the generator of the matrix zonotope
    Z_K = cell(Q{1}.gens,1);
    
    for j = 1:Q{1}.gens
        for i = 1:length(Q)
            Q_{i} = Q{i}.generator{j}; 
        end
        
        temp = quadMap(Z,Q_);
        Z_K{j} = [temp.c,temp.G];
    end
    
    % overall zonotope
    Zquad = zonotope(Z_D.c, [Z_D.G,[Z_K{:}]]);
    Zquad = compact_(Zquad,'zeros',eps);
end

% ------------------------------ END OF CODE ------------------------------
