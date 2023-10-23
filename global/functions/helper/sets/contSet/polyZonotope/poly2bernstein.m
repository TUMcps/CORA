function [B,E] = poly2bernstein(G,E,dom)
% poly2bernstein - Convert a polynomial to a bernstein polynomial
%
% Syntax:
%    [B,E] = poly2bernstein(G,E,dom)
%
% Inputs:
%    G - generator matrix containing the coefficients of the polynomial
%    E - exponent matrix containing the exponents of the polynomial
%    dom - interval that represents the domain for each variable
%
% Outputs:
%    B - coefficients of the bernstein polynomial
%    E - exponent matrix of the bernstein polnomial
%
% Example: 
%   E = [1 0 2;0 1 1];
%   G = [1 2 3];
%   dom = interval([-1;-1],[1;1]);
%    
%   [B,E_] = poly2bernstein(G,E,dom);
%
% References: 
%   [1] S. Ray et al. "A Matrix Method for Efficient Computation of
%       Bernstein Coefficients"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm/interval, polyZonotope/interval

% Authors:       Niklas Kochdumper
% Written:       03-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    % Implementation of the "Bernstein_Matrix" algorithm from [1]

    l = size(E,1);
    n = max(max(E));
    h = size(G,1);
    
    % Preprocessing: compute coefficient matrix A
    len = (n+1)^(l-1);
    A = repmat({zeros(n+1,len)},[h,1]);
    
    for i = 1:size(E,2)
        ind = aux_getIndices(E(:,i),len,n);
        for j = 1:h
            A{j}(ind(1),ind(2)) = G(j,i);
        end
    end

    % Step 1: Compute the binomial coefficients
    C = aux_Binomial_coefficient(n);
    
    % Step 2: Compute inverse of U_x
    Ux = aux_InverseUx(n,C);
    
    % Step 3: Iterate
    M = cell(l,1);
    infi = infimum(dom);
    sup = supremum(dom);
    
    for r = 1:l 
        if r == 1 || infi(r) ~= infi(1) || sup(r) ~= sup(1) 
            
            % Compute inverse of V_x
            Vx = aux_InverseVx(n,dom(r));

            % Compute inverse of W_x
            Wx = aux_InverseWx(n,dom(r),C);

            % Product of all the inverse matrices
            M{r} = Ux * Vx * Wx;
        
        else
            M{r} = M{1}; 
        end
    end
    
    % Step 4: Iterate
    for j = 1:h
        for r = 1:l       
            A{j} = aux_transposeMatrix(A{j},r,len,l,n);
            A{j} = M{r}*A{j};
            A{j} = aux_transposeMatrix(A{j},r,len,l,n);
        end
    end
    
    % assemble output arguments
    if nargout > 1
        
        % construct exponent matrix for bernstein polynomials
        E = combinator(n+1,l)' - 1;
        B = zeros(h,size(E,2));
        counter = 1;
        
        % construct vector of bernstein coefficients
        for i = 1:length(B)
            ind = aux_getIndices(E(:,i),len,n);
            for j = 1:h
                B(j,counter) = A{j}(ind(1),ind(2));
            end
            counter = counter + 1;
        end
    else
        if h == 1
            B = A{1};
        else
            B = A; 
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function C = aux_Binomial_coefficient(n)
% Implementation of the "Binomial_coefficient" algorithm in [1]

    C = NaN * ones(n+1,n+1);
    C(2,1) = 1;
    
    for i = 2:n
        
       C(i,i) = 1;
       C(i+1,1) = 1;
       
       for k = 1:(i-1)
          
           C(i+1,k+1) = i/(i-k)*C(i,k+1);
           
       end
    end
end

function Ux = aux_InverseUx(n,C)
% Implementation of the "InverseUx" algorithm in [1]

    Ux = zeros(n+1,n+1);
    
    Ux(1:end,1) = ones(n+1,1);
    Ux(n+1,2:end) = ones(1,n);
    
    for i = 1:n-1
       
        for j = 1:i
            
            Ux(i+1,j+1) = C(i+1,i-j+1)/C(n+1,j+1);
   
        end
    end
end

function Vx = aux_InverseVx(n,x)
% Implementation of the "InverseVx" algorithm in [1]

    wid_x = 2*rad(x);
    
    Vx = zeros(n+1,n+1);
    Vx(1,1) = 1;
    
    wid_x_ = zeros(n+1,1);
    wid_x_(1) = 1;
    
    for i = 1:n
        wid_x_(i+1) = wid_x_(i)*wid_x;
        Vx(i+1,i+1) = wid_x_(i+1);
    end
end

function Wx = aux_InverseWx(n,x,C)
% Implementation of the "InverseWx" algorithm in [1]

    % All the powers of the infimum
    inf_x = infimum(x);
    
    inf_x_ = zeros(n+1,1);
    inf_x_(1) = 1;
    
    for i = 1:n
       inf_x_(i+1) = inf_x_(i) * inf_x;
    end

    % Construction of the inverse matrix Wx
    Wx = zeros(n+1,n+1);
    
    Wx(1,1) = 1;
    
    for i = 0:n-1
        
        for j = i+1:n
           Wx(i+1,j+1) = C(j+1,i+1) * inf_x_(j-i+1);
        end
        
        Wx(i+2,i+2) = 1;
    end 
end

function ind = aux_getIndices(e,len,n)
% get the matrix index for the exponent vector e

    % column index = exponent of first variable
    ind = zeros(2,1);
    ind(1) = e(1) + 1;
    
    % row index -> loop over all remaining variables 
    for i = length(e):-1:2
        len_ = len/(n+1);
        ind(2) = ind(2) + e(i) * len_;
        len = len_;
    end
    
    ind(2) = ind(2) + 1;
end

function A = aux_transposeMatrix(A,r,len,l,n)
% Implementation of matrix transposition as required for the evaluation of
% equation (63) in [1]

    if r == 2
        A = reshape(permute(reshape(A,[n+1,n+1,len/(n+1)]),[2,1,3]),[n+1,len]);
    elseif r ~= 1
        if r == l
            A = aux_transpose_(A);
        else
            len_ = (n+1)^(r-1);
            temp = mat2cell(A,n+1,len_ * ones(1,len/len_));
            temp = cellfun(@aux_transpose_,temp,'UniformOutput',false);
            A = [temp{:}];
        end
    end
end

function A = aux_transpose_(A)

    n = size(A,1);
    p = size(A,2)/n;
    
    temp = reshape(permute(reshape(A,[n,p,n]),[2,1,3]),[1,n*p,n]);
    A = permute(temp,[3,2,1]);
end

% ------------------------------ END OF CODE ------------------------------
