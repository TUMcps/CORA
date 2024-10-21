function SpS = generateESumRep(SpS)
% generateESumRep - generates the Existential Sum Representation of a
%    spectrahedral shadow and saves it in the hidden property "ESumRep"
%
% Syntax:
%    SpS = generateESumRep(SpS)
%
% Inputs:
%    SpS - spectraShadow object with empty field "ESumRep"
%
% Outputs:
%    SpS - spectraShadow object with non-empty field "ESumRep"
% 
% Example:
%    SpS = spectraShadow([eye(3) eye(3) eye(3)],1,[3 1]);
%    generateESumRep(SpS);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       05-May-2023 
% Last update:   02-August-2023 (AK, implemented degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isempty(SpS.ESumRep.val{1}) || ~isempty(SpS.ESumRep.val{2})
    return;
end

% get projection parameters
c = SpS.c;


% remove translation should it exist
if ~all(withinTol(c,zeros(size(c)),1e-12))
    SpS_trans = aux_elimTrans(SpS);
else
    SpS_trans = SpS;
end

G = SpS_trans.G;

% We prepare the ground by computing the svd of G:
[U,Sigma,V] = svd(full(G));
r = rank(full(G));
% Isolating the non-zero part of Sigma
D = Sigma(1:r,1:r);
n = size(Sigma,1);
m = size(Sigma,2);


% as well as the 'inner' spectraShadow, i.e., the one that is not
% projected:
SpS_inner = spectraShadow(SpS_trans.A);
SpS_inner.ESumRep.val = {SpS_trans.A []};
% We want to transform S_inner, such that T(S_inner) = S, but in such a way
% that we keep track of ESumRep
% So first, we need to rotate S_inner by V':
SpS_inner = priv_mtimesIsomorphism(SpS_inner, V');
[A0, Ai] = priv_getCoeffMatrices(spectraShadow(SpS_inner.ESumRep.val{1}));

% Separating S_inner.A into an ESumRep that decides which variables
% contribute, and which will be 'discounted' by the zero-part in S
ESumRep{1} = [A0 cat(2,Ai{1:r})];
ESumRep{2} = cat(2,Ai{r+1:end});
% Define a temporary spectrahedron, for the 'upper' part affected by the
% non-zero part of Sigma
SpS_temp = spectraShadow(ESumRep{1});
SpS_temp.ESumRep.val = ESumRep;
% Now multiply it by D
SpS_temp = priv_mtimesIsomorphism(SpS_temp, D);
SpS_inner = aux_cartProdZeros(SpS_temp,n-r);

% And finally, we just multiply S_inner again with U. This should give us
% back our S, but this time with the existential sum representation
% computed along the way
SpS_existential = priv_mtimesIsomorphism(SpS_inner, U);
SpS.ESumRep.val = SpS_existential.ESumRep.val;

end


% Auxiliary functions -----------------------------------------------------

function SpS_adj = aux_elimTrans(SpS_t)
% eliminate translation parameter of a spectrahedron by implementing it
% into the coefficient matrices    

    G = SpS_t.G;
    c = SpS_t.c;
    
    m = size(G,2);
    
    % Temporary fix to prevent numerical errors; if someone finds a better
    % solution, be free to change this :)
    G(abs(G)<= 1e-5) = 0;
    
    % There are two methods how one can eliminate the need for a
    % translation vector. Which one is better depends on whether G can
    % represent c
    if rank(full(G)) == rank([full(G) full(c)])
        % So, if G is not necessarily full rank, but at least can represent
        % c, we can find some coordinates x such that Gx = c
        x = G \ c;
        
        % And so we can continue with 'eliminating' x
    
    
        % get separated coefficient matrices
        [A0, Ai] = priv_getCoeffMatrices(SpS_t);

        % compute new first coefficient 
        for i = 1:m
            A0 = A0 - x(i)*Ai{i};
        end

        % set up new A
        A = [A0 cat(2,Ai{:})];

        % generate new spectraShadow object
        SpS_adj = spectraShadow(A,zeros(size(c)),G);
    else
        % If that is not the case, we are facing a degenerate
        % spectrahedron, and things are a bit more complicated.
        % Basically, we need to represent S as
        % S = {[G c]x_tilde | ...},
        % where x_tilde = [x x_end]^T, and x_end is constrained in such a
        % way that only x_end = 1 is allowed.
        
        % get separated coefficient matrices
        [A0, Ai] = priv_getCoeffMatrices(SpS_t);
        
        A0_new = blkdiag(A0,sparse([1 0;0 -1]));
        Ai_new = cell([1 size(G,2)+1]);
        for i = 1:size(G,2)
            Ai_new{i} = blkdiag(Ai{i},sparse(2,2));
        end
        Ai_new{end} = blkdiag(sparse(size(A0,1), size(A0,2)), sparse([-1 0;0 1]));
        
        A = [A0_new cat(2,Ai_new{:})];
        
        SpS_adj = spectraShadow(A,zeros(size(c)),[G c]);
    end
end

function SpS_0 = aux_cartProdZeros(SpS,d)
% Computes the Cartesian product S x {0}^d, but only by adding constraints
% to ESumRep
    c = SpS.c;
    G = SpS.G;
    A = SpS.A;
    ESumRep = SpS.ESumRep.val;
    
    c_0 = [c;zeros([d 1])];
    G_0 = [G;sparse(d,size(G,2))];
    A_0 = A;
    
    % Extend ESumRep{1}
    k = size(ESumRep{1},1);
    m = size(ESumRep{1},2) / k - 1;
    
    % We first extend the existing matrices
    ESumRep_0_B0 = blkdiag(ESumRep{1}(:,1:k), sparse(2*d,2*d));
    ESumRep_0_Bi = cell([1 m]);
    for i=1:m
        ESumRep_0_Bi{i} = blkdiag(ESumRep{1}(:,(i*k + 1):(i*k+k)),sparse(2*d,2*d));
    end
    ESumRep_0{1} = [ESumRep_0_B0 cat(2,ESumRep_0_Bi{:})];
    
    % Now we add the constraints that also add new variables, which are all
    % supposed to be zero
    ESumRep_0_Bi_additional = cell([1 d]);
    for i=1:d
        v = zeros([k+2*d 1]);
        v(k+2*(i-1)+1) = 1;
        v(k+2*i) = -1;
        ESumRep_0_Bi_additional{i} = spdiags(v,0,k+2*d,k+2*d);
    end
    ESumRep_0{1} = [ESumRep_0{1} cat(2,ESumRep_0_Bi_additional{:})];
    
    % Extend ESumRep{2}. We only need to deal with the 'padding'
    k = size(ESumRep{2},1);
    m = size(ESumRep{2},2) / k;
    
    if k==0
        ESumRep_0{2} = [];
    else
        ESumRep_0_Ci = cell([1 m]);
        for i=1:m
            ESumRep_0_Ci{i} = blkdiag(ESumRep{2}(:,((i-1)*k + 1):((i-1)*k+k)),sparse(2*d,2*d));
        end
        ESumRep_0{2} = cat(2,ESumRep_0_Ci{:});
    end
    
    SpS_0 = spectraShadow(A_0, c_0, G_0);
    SpS_0.ESumRep.val = ESumRep_0;
end

% ------------------------------ END OF CODE ------------------------------
