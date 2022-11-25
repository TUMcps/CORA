function res = testLongDuration_bernsteinPoly
% testLongDuration_bernsteinPoly - unit-test that tests for the conversion
%    of a polynomial to a bernstein polynomial as implemented by the
%    function "poly2bernstein"
%
% Syntax:  
%    res = testLongDuration_bernsteinPoly
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: poly2bernstein

% Author:       Niklas Kochdumper
% Written:      31-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    N = 5;
    M = 10;

    % loop over all dimensions
    for n = 1:N
        
       % loop over all test-cases
       for i = 1:M
          
           % generate random polynomial zonotope
           pZ = polyZonotope.generateRandom(1,[],n);
           
           n_ = size(pZ.expMat,1);
           expMat = [zeros(n_,1),pZ.expMat];
           G = [pZ.c,pZ.G];
           
           dom = interval(-ones(n_,1),ones(n_,1));
           
           % convert to bernstein polynomial
           [B,expMat_] = poly2bernstein(G,expMat,dom);
           
           % test for random points if the two polynomials are identical
           for j = 1:M
               
              p = randPoint(dom);
              
              val = regularPolynom(p,G,expMat);
              val_ = bernsteinPolynom(p,B,expMat_);
              
              if abs(val-val_) > 1e-13
                 error('Unit-test "testLongDuration_bernsteinPoly" failed!'); 
              end
           end
       end
    end
    
    res = 1;
end


% Auxiliary Functions -----------------------------------------------------

function val = regularPolynom(x,G,expMat)
% evaluate a polynomial at the point x

    val = 0;
    
    for i = 1:length(G)
       temp = 1;
       for j = 1:size(expMat,1)
          temp = temp * x(j)^expMat(j,i); 
       end
       val = val + G(i) * temp;
    end
end

function val = bernsteinPolynom(x,B,expMat)
% evaluate a bernstein polynomial at point x

    l = length(x);
    n = max(max(expMat));
    
    dom = interval(-ones(l,1),ones(l,1));
    
    val = 0;
    
    for i = 1:length(B)
       val = val + B(i) * bernsteinPolynom_(x,dom,n,expMat(:,i)); 
    end
end

function val = bernsteinPolynom_(x,dom,n,e)
% evaluate a single bernstein-basis polynomial at point x

    infi = infimum(dom);
    sup = supremum(dom);
    
    val = 1;
    
    for i = 1:length(x)
       
        val = val * nchoosek(n,e(i)) * (x(i)-infi(i))^e(i) * ...
              (sup(i)-x(i))^(n-e(i))/(sup(i)-infi(i))^n;
    end
end

%------------- END OF CODE --------------