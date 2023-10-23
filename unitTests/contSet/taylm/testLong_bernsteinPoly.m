function res = testLong_bernsteinPoly
% testLong_bernsteinPoly - unit-test that tests for the conversion
%    of a polynomial to a bernstein polynomial as implemented by the
%    function "poly2bernstein"
%
% Syntax:
%    res = testLong_bernsteinPoly
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: poly2bernstein

% Authors:       Niklas Kochdumper
% Written:       31-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % assume true
    res = true;

    N = 5;
    M = 10;

    % loop over all dimensions
    for n = 1:N
        
       % loop over all test-cases
       for i = 1:M
          
           % generate random polynomial zonotope
           pZ = polyZonotope.generateRandom('Dimension',1,'NrFactors',n);
           
           n_ = size(pZ.E,1);
           E = [zeros(n_,1),pZ.E];
           G = [pZ.c,pZ.G];
           
           dom = interval(-ones(n_,1),ones(n_,1));
           
           % convert to bernstein polynomial
           [B,E_] = poly2bernstein(G,E,dom);
           
           % test for random points if the two polynomials are identical
           for j = 1:M
               
              p = randPoint(dom);
              
              val = aux_regularPolynom(p,G,E);
              val_ = aux_bernsteinPolynom(p,B,E_);
              
              if ~withinTol(val,val_,1e-13)
                 res = false; return
              end
           end
       end
    end
    
end


% Auxiliary functions -----------------------------------------------------

function val = aux_regularPolynom(x,G,E)
% evaluate a polynomial at the point x

    val = 0;
    
    for i = 1:length(G)
       temp = 1;
       for j = 1:size(E,1)
          temp = temp * x(j)^E(j,i); 
       end
       val = val + G(i) * temp;
    end
end

function val = aux_bernsteinPolynom(x,B,E)
% evaluate a bernstein polynomial at point x

    l = length(x);
    n = max(max(E));
    
    dom = interval(-ones(l,1),ones(l,1));
    
    val = 0;
    
    for i = 1:length(B)
       val = val + B(i) * aux_bernsteinPolynom_(x,dom,n,E(:,i)); 
    end
end

function val = aux_bernsteinPolynom_(x,dom,n,e)
% evaluate a single bernstein-basis polynomial at point x

    lb = infimum(dom);
    ub = supremum(dom);
    
    val = 1;
    
    for i = 1:length(x)
       
        val = val * nchoosek(n,e(i)) * (x(i)-lb(i))^e(i) * ...
              (ub(i)-x(i))^(n-e(i))/(ub(i)-lb(i))^n;
    end
end

% ------------------------------ END OF CODE ------------------------------
