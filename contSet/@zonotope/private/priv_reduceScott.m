function Zred = priv_reduceScott(Z,order)
% priv_reduceScott - reduce zonotope so that its order stays below a
%    specified limit according to [1]. This reduction method is especially
%    suited for constrained zonotopes
%
% Syntax:
%    Zred = priv_reduceScott(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/reduce
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper
% Written:       16-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine the number of zonotope generators that get reduced
n = dim(Z);
N = ceil(size(Z.G,2) - order * n);

% check if it is necessary to reduce the order
if N <= 0
    Zred = Z;
    return
end

% compute low echelon form of the generator matrix using 
% Gauss-Jordan elimination
[A,indPer,fullRank] = aux_rrefInfty(Z.G);

% if matrix has rank deficit, reduce with combastel as a back-up strategy
if ~fullRank        
   Zred = reduce(Z,'combastel',order);
   return
end

% ...matrix is full rank      

% reorder columns of generator matrix to the form [T V] 
% according to the low echelon form [I R]
T = Z.G(:,indPer(1:n));

R = A;
R(:,1:n) = [];

% loop over all generators that get removed
for i = 1:N

    % compute costs (= volume error) for each generator in V
    costs = zeros(size(R,2),1);

    for j = 1:length(costs)
        r = R(:,j);
        vol = 1 + sum(abs(r));
        vol_ = prod(abs(r) + ones(n,1));
        costs(j) = vol_ - vol;              
    end

    % determine generator with minimal costs
    [~,ind] = min(costs);
    r = R(:,ind);

    % update generator matrices R and T
    T = T * (eye(n) + diag(abs(r)));

    R_ = R;
    R_(:,ind) = [];
    R = diag(1./(ones(n,1)+abs(r))) * R_;            
end

% recover the resulting reduced zonotope
Zred = zonotope(Z.c, [T, T*R]);

end  
    

% Auxiliary functions -----------------------------------------------------
    
function [A,indPer,fullRank] = aux_rrefInfty(A)
% Transform the matrix to Reduced Echelon Form using Gauss-Jordan
% elimination with full pivoting. The row elements with the largest 
% absolute values relative to the infinity norm of their row are chosen as 
% the pivot elements.

    [m,n] = size(A);
    fullRank = true;
    
    % compute the default tolerance
    tol = max(m,n)*eps(class(A))*norm(A,'inf');

    % loop over the entire matrix
    indPer = 1:n;

    for i = 1:m
        
       % unreduced submatrix
       Atemp = A(i:m,i:n);
        
       % divide each row by it's inifinity norm
       normInf = sum(abs(Atemp),2);
       Atemp = diag(1./normInf) * Atemp;
        
       % find value and index of largest element in the remainder of column j
       if size(Atemp,1) > 1
           [valTemp,indTemp] = max(abs(Atemp));
           [p,indC] = max(valTemp);
           indR = indTemp(indC) + i - 1;
           indC = indC + i - 1;
       else
          [p,indC] = max(abs(Atemp));
          indC = indC + i - 1;
          indR = i;
       end
       
       if p <= tol || isnan(p)
           fullRank = false;
       else
          
          % bring the row with the pivot element up
          A([i indR],:) = A([indR i],:);
          
          % bring the column with the pivot element to the front
          A(:,[i indC]) = A(:,[indC i]);
          temp = indPer(i);
          indPer(i) = indPer(indC);
          indPer(indC) = temp;
          
          % divide the pivot row by the pivot element
          Ai = A(i,:)/A(i,i);  
          
          % subtract multiples of the pivot row from all the other rows
          A(:,i:n) = A(:,i:n) - A(:,i)*Ai(:,i:n);
          A(i,i:n) = Ai(i:n);
       end
    end
end

% ------------------------------ END OF CODE ------------------------------
