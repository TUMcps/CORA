function res = rescale(obj,varargin)
% resclae - Rescales the domains for the factors ksi of a constrained
%           zonotope
%
% Syntax:  
%    res = rescale(obj)
%    res = rescale(obj, method)
%
% Inputs:
%    obj - c-zonotope object
%    method - method used to determine the tighend domain for the zonotope
%             factors ('exact' or 'iter')
%
% Outputs:
%    res - c-zonotope object
%
% Example:
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    cRes = rescale(cZono,'iter');
%    
%    hold on
%    plotZono(cZono,[1,2]);
%    plotZono(cRes,[1,2],{'g'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      18-December-2017
% Last update:  12-May-2018
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    method = 'exact';
    if nargin == 2
       method = varargin{1};
    end

    % determine the new domain for each factor ksi
    if isempty(obj.A)
        error('No constrains left!');
    elseif isempty(obj)
        error('Set is empty!');
    elseif isempty(obj.ksi)
        if strcmp(method,'iter')
            obj = preconditioning(obj);
            [domKsi, R] = ksi_iterative( obj );
        else
            domKsi = ksi_optimizer( obj );
        end

        ksi_l = infimum(domKsi);
        ksi_u = supremum(domKsi);
    else
        ksi_l = min(obj.ksi,[],2);
        ksi_u = max(obj.ksi,[],2);
    end

    % new basis
    ksi_m = (ksi_u + ksi_l)/2;
    ksi_r = (ksi_u - ksi_l)/2;

    % rescale R (Equation (A.3) in reference paper [1])
    if strcmp(method,'iter')
        rho_l = infimum(R);
        rho_u = supremum(R);

        obj.R = [(rho_l - ksi_m)./ksi_r, (rho_u - ksi_m)./ksi_r];
    end

    % rescale c-zonotope (Equation (24) in reference paper [1])
    temp = diag(ksi_r);

    G = obj.Z(:, 2:end);
    c = obj.Z(:, 1) + G*ksi_m;
    obj.Z = [c, G * temp];
    obj.b = obj.b - obj.A * ksi_m;
    obj.A = obj.A * temp;

    if ~isempty(obj.ksi)
       temp = ones(1,size(obj.ksi,2));
       ksi = obj.ksi - ksi_m*temp;
       obj.ksi = zeros(size(ksi)) + ksi.*((1./ksi_r)*temp);
    end

    % output the result
    res = obj;

end




% Auxiliary Functions -----------------------------------------------------

function obj = preconditioning(obj)
% preconditioning of the constraint matrix to improve rescaling    

    % bring constraint matrices to Reduced Echelon Form
    [obj.A,obj.b,indPer] = rrefInfty(obj.A,obj.b);
    
    % adapt the generator matrix to the new order of factors
    c = obj.Z(:,1);
    G = obj.Z(:,2:end);
    
    obj.Z = [c, G(:,indPer)];
end


function [A_,b_,indPer] = rrefInfty(A,b)
% Transform the matrix [A b] to Reduced Echelon Form using Gauss-Jordan
% elimination with full pivoting. The row elements with the largest 
% absolute values relative to the inifinity norm of their row are chosen as 
% the pivot elements.

    A_ = [A,b];
    [m,n] = size(A_);
    
    % compute the default tolerance
    tol = max(m,n)*eps(class(A_))*norm(A_,'inf');

    % loop over the entire matrix
    indPer = 1:n-1;

    for i = 1:m
        
       % unreduced submatrix
       Atemp = A_(i:m,i:n-1);
        
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
       
       if (p > tol)
          
          % bring the row with the pivot element up
          A_([i indR],:) = A_([indR i],:);
          
          % bring the column with the pivot element to the front
          A_(:,[i indC]) = A_(:,[indC i]);
          temp = indPer(i);
          indPer(i) = indPer(indC);
          indPer(indC) = temp;
          
          % divide the pivot row by the pivot element
          Ai = A_(i,:)/A_(i,i);  
          
          % subtract multiples of the pivot row from all the other rows
          A_(:,i:n) = A_(:,i:n) - A_(:,i)*Ai(:,i:n);
          A_(i,i:n) = Ai(i:n);
       end
    end
    
    % assign output arguments
    b_ = A_(:,end);
    A_ = A_(:,1:n-1);
end

%------------- END OF CODE --------------