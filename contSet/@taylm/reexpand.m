function res = reexpand(obj,int)
% reexpand - Reexpand a taylor model at the center of the domain int
%
% Syntax:  
%    res = reexpand(obj, int)
%
% Inputs:
%    obj - taylm object
%    int - new domain of the taylor model (class interval)
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval
% Subfunctions: intPower, intMul, evalInt
% MAT-files required: none
%
% See also: taylm

% Author:       Niklas Kochdumper
% Written:      14-April-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE -------------

    % calculate the new expansion point
    x_e = mid(int);
    
    % calculate the scale factor between the old and the new interval for
    % all dimensions
    scaleFact = rad(int);
    
    % caclulate the taylor series of the taylor model around the new
    % expansion point
    K = max(max(obj.monomials));
    monAll = monomials(length(x_e),K);
    
    counter = 1;
    coeff = zeros(1,size(monAll,1));
    monNew = zeros(size(monAll));
    
    for i = 1:size(monAll,1)
        
        % calculate derivative of the talyor model
        e = monAll(i,:);
        coeff(counter) = evalDiff(obj,x_e,e);
        
        if coeff(counter) ~= 0
            % scale the domain of the variables to [-1,1]
            coeff(counter) = coeff(counter) * prod(scaleFact.^(e'));
            monNew(counter,:) = e;
            counter = counter + 1;
        end
    end
    
    coeff = coeff(1:counter-1);
    monNew = monNew(1:counter-1,:);

    % construct the new taylor model object
    res = obj;
    res.coefficients = coeff';
    res.monomials = [sum(monNew,2),monNew];
end



function c = evalDiff(obj,x,diff)
% evaluate the derivative of a taylor model

    c = 0;
    
    for i = 1:length(obj.coefficients)
       e = obj.monomials(i,2:end);
       e_ = e-diff;
       if e_ >= 0
           cTemp = 1;
           for j = 1:length(e)
               cTemp = cTemp * prod(e_(j) + 1 : e(j));
           end

           c = c + obj.coefficients(i) * cTemp * prod(x.^(e_'));
       end
    end

    c = c * 1/prod(factorial(diff));
end

function mon = monomials(M,N)
% construct all monomials up to the degree N

    % inititial monomial
    mon = zeros(1,M);
    
    % loop over all monomials up to a degree of N
    while sum(mon(end,:)) <= N
        
       temp = nextMonomial(M,mon(end,:)); 
       mon = [mon;temp];     
    end
    
    mon = mon(1:end-1,:);

end

function x = nextMonomial(M,x)
% calculte the next monomial according to the graded lexicographic order

      j = 1;

      for i = 2 : M
        if ( 0 < x(i) )
          j = i;
          break
        end
      end

      if ( j == 1 )
        t = x(1);
        x(1) = 0;
        x(M) = t + 1;
      elseif ( j < M )
        x(j) = x(j) - 1;
        t = x(1) + 1;
        x(1) = 0;
        x(j-1) = x(j-1) + t;
      elseif ( j == M )
        t = x(1);
        x(1) = 0;
        x(j-1) = t + 1;
        x(j) = x(j) - 1;
      end
end

%------------- END OF CODE --------------
