function Zred = reduceScott(Z,order)
% reduceScott - Reduce zonotope so that its order stays below a specified
%               limit. This reduction method is especially suited for
%               constrained zonotopes
%
% Syntax:  
%    Zred = reduceScott(Z,order)
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

% Author:       Niklas Kochdumper
% Written:      16-July-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % Implementation of the zonotope reduction technique desribed in the
    % Appendix of [1]
    
    % extract center and generators
    c = Z.Z(:,1);
    G = Z.Z(:,2:end);
    
    % determine the number of zonotope generators that get reduced
    n = size(G,1);
    N = ceil(size(G,2) - order * n);

    % check if it is necessary to reduce the order
    if N > 0
    
        % compute low echelon form of the generator matrix using 
        % Gauss-Jordan elimination
        [A,jb] = rrefInfty(G);

        % check matrix rank
        if jb ~= size(G,1)         % matrix has rank defincit

           % reduce with combastel as a back-up strategy
           Zred = reduce(Z,'combastel',order);


        else                       % matrix is full rank

            % reorder columns of generator matrix to the form [T V] 
            % according to the low echelon form [I R]
            T = G(:,jb);

            R = A;
            R(:,jb) = [];

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
            Zred = zonotope([c, T, T*R]); 
        end
    
    else
        Zred = Z;
    end
end  
    

% Auxiliary Functions -----------------------------------------------------
    
function [A,jb] = rrefInfty(A)
% Transform a matrix to Reduced Echelon Form using Gauss-Jordan
% elimination. The row elements with the largest absolute values are chosen
% as the pivot elements.

    [m,n] = size(A);
    
    % compute the default tolerance
    tol = max(m,n)*eps(class(A))*norm(A,'inf');
    
    % divide each row by it's inifinity norm
    normInf = sum(abs(A),2);
    A = diag(normInf) * A;

    % loop over the entire matrix
    i = 1;
    j = 1;
    jb = [];

    while (i <= m) && (j <= n)
        
       % find value and index of largest element in the remainder of column j
       [p,k] = max(abs(A(i:m,j))); k = k+i-1;
       
       if (p <= tol)
           
          % the column is negligible, zero it out
          A(i:m,j) = 0; 
          j = j + 1;
          
       else
           
          % remember column index
          jb = [jb j];
          
          % swap i-th and k-th rows
          A([i k],j:n) = A([k i],j:n);
          
          % divide the pivot row by the pivot element
          Ai = A(i,j:n)/A(i,j);  
          
          % subtract multiples of the pivot row from all the other rows
          A(:,j:n) = A(:,j:n) - A(:,j)*Ai;
          A(i,j:n) = Ai;
          i = i + 1;
          j = j + 1;
       end
    end
end

%------------- END OF CODE --------------