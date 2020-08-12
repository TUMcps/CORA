function res = test_zonotope_polytope
% test_polytope - unit test function of polytope
%
% Syntax:  
%    res = test_zonotope_polytope
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Analytical Test ---------------------------------------------------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain polytope
P = polytope(Z1);

% extract MPT object
P = get(P,'P');

% obtain halfspace matrix
C = P.H(:,1:end-1);

% obtain distance vector
d = P.H(:,end);

% true results
true_C = [ 0.554700196225229   0.832050294337844; ...
           0.832050294337844   0.554700196225229; ...
           0.970142500145332   0.242535625036333; ...
          -0.554700196225229  -0.832050294337844; ...
          -0.832050294337844  -0.554700196225229; ...
          -0.970142500145332  -0.242535625036333];
      
true_d =  [2.773500981126146; 0; 0; 5.547001962252290; 5.547001962252290; 7.276068751089989];

% check result
for i = size(C,1)
    [temp,ind] = ismembertol(C(i,:),true_C,1e-13,'ByRows',true);
    if ~temp || abs(d(i)-true_d(ind)) > 1e-13 || size(C,1) ~= size(true_C,1)
        error('Analytical test failed!'); 
    end
end



% Random Test -------------------------------------------------------------

res = 1;
N = 100;
Tol = 1e-12;

% loop over different dimensions
for n = 1:4
    
   % loop over the number of generators
   for m = 1:6
      
       % construct random zonotope
       zono = zonotope(rand(n,m+1)-0.5*ones(n,m+1));
       
       % convert to polytope and extract inequality constraints
       poly = polytope(zono);
       P = get(poly,'P');
       A = P.A;
       b = P.b;
       
       % draw random points from inside the zonotope
       pointsIn = zeros(n,N);
       
       for i = 1:N
          pointsIn(:,i) = randPoint(zono); 
       end
       
       % draw random points outside the zonotope
       pointsOut = zeros(n,N);
       inter = interval(zono);
       inter_ = interval(center(inter)-2*rad(inter),center(inter)+2*rad(inter));
       zono_ = zonotope(inter_);
       counter = 1;
       
       while counter <= N
           p = randPoint(zono_);
           if ~containsPoint(inter,p)
              pointsOut(:,counter) = p;
              counter = counter + 1;
           end
       end
       
       % check if the inside points fullfil the inequality constraints
       temp = A*pointsIn - b * ones(1,N);
       
       if any(any(temp > Tol))
          res = 0;
          error('Random test failed! '); 
       end
       
       % check if the outside points violate the inequality constraints
       temp = A*pointsOut - b * ones(1,N);
       
       for i = 1:N
           if ~any(temp(:,i) > -Tol)
              res = 0;
              error('Random test failed!'); 
           end  
       end
   end
end




if res
    disp('test_zonotope_polytope successful');
end

%------------- END OF CODE --------------
