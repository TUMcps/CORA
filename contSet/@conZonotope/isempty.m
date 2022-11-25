function res = isempty(obj)
% isempty - returns 1 if a constrained zonotope is empty and 0 otherwise
%
% Syntax:  
%    res = isempty(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%   res - result in {0,1}
%
% Example:
%    c = [0;0];
%    G = [1 0 1;0 1 1];
%    A = [1 1 1];
%    b = 4;
%    cZ = conZonotope([c,G],A,b);
%
%    isempty(cZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      15-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check if the zonotope is empty
res = isempty(zonotope(obj.Z));

% check if the constraints are satisfiable 
if ~res && ~isempty(obj.A)
   
   % Calculate null space of the constraints
   Neq=null(obj.A);   
   
   % Calculate a single point that satisfies the constraints
   x0=pinv(obj.A)*obj.b;
   
   % Define tolerance
   Tol = 1e-10;

   if norm(obj.A*x0-obj.b)>Tol*norm(obj.b)  % infeasible

      res = 1;

   elseif isempty(Neq)  % null space empty -> set is a single point

       % construct the inequatility constraints (unit cube)
       n = size(obj.Z,2)-1;
       A = [eye(n);-eye(n)];
       b = [ones(n,1);ones(n,1)];
     
       % check if the point satisfies the inequality constraints 
       if ~all(A*x0<=b)

          res = 1;
       end
       
   else     % check if the null-space intersects the unit-cube
      
       temp = ones(size(obj.A,2),1);
       unitCube = interval(-temp,temp);
       
       % loop over all constraints (= hyperplanes)
       for i = 1:size(obj.A,1)
           
          % hyperplane from a constraint does not intersect the unit cube
          % -> set is empty
          if ~isIntersecting(halfspace(obj.A(i,:),obj.b(i)),unitCube)
              res = 1;
              return ;
          end
       end
       
       % use linear programming to check if the constrained zonotope is
       % empty (this seems to be more robust than the previous solution
       % using the mptPolytope/isempty function)
       if size(obj.A,1) >= 1
           
           options = optimoptions('linprog','display','off', ...
                               'OptimalityTolerance',1e-10);	
                           
           p = size(obj.A,2);
           ub = ones(p,1); lb = -ub; f = ones(p,1);
           
           [~,~,exitflag] = linprog(f,[],[],obj.A,obj.b,lb,ub,options);

            res = 0;

            if exitflag == -2
                res = 1; 
            end
       end
              
   end 
end

%------------- END OF CODE --------------