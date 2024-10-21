function res = testLong_zonotope_polytope
% testLong_zonotope_polytope - unit test function of conversion to polytope
%
% Syntax:
%    res = testLong_zonotope_polytope
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
% See also: none

% Authors:       Matthias Althoff
% Written:       26-July-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

N = 100;
Tol = 1e-12;

% loop over different dimensions
for n = 1:4
    
    % loop over the number of generators
    for m = 1:6
        
        % construct random zonotope
        Z = zonotope(rand(n,m+1)-0.5*ones(n,m+1));
        
        % convert to polytope and extract inequality constraints
        poly = polytope(Z);
        A = poly.A;
        b = poly.b;
        
        % draw random points from inside the zonotope
        pointsIn = zeros(n,N);
        
        for i = 1:N
            pointsIn(:,i) = randPoint(Z); 
        end
        
        % draw random points outside the zonotope
        pointsOut = zeros(n,N);
        I = interval(Z);
        I_ = interval(center(I)-2*rad(I),center(I)+2*rad(I));
        Z_ = zonotope(I_);
        counter = 1;
        
        while counter <= N
            p = randPoint(Z_);
            if ~contains(I,p)
                pointsOut(:,counter) = p;
                counter = counter + 1;
            end
        end
        
        % check if the inside points fulfill the inequality constraints
        margin = A*pointsIn - b * ones(1,N);
        
        assertLoop(margin <= 0 | withinTol(margin,0,Tol),n,m);
        
        % check if the outside points violate the inequality constraints
        margin = A*pointsOut - b * ones(1,N);
       
        for i = 1:N
            assertLoop(any(margin(:,i) > 0 | withinTol(margin(:,i),0,Tol)),m,n,i);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
