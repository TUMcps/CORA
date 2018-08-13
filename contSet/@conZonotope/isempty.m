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
%    ---
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

% check if the zonotope is empty by calling the superclass method
res = isempty@zonotope(obj);

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
       
       % check if the null-space intersects the unit-cube
       if size(obj.A,1) >= 1
           
           % construct inequality constraints for the unit cube
           n = size(obj.Z,2)-1;
           A = [eye(n);-eye(n)];
           b = [ones(n,1);ones(n,1)];

           % transform the constraints to the null space
           A_ = A*Neq;
           b_ = b-A*x0;

           % create mptPolytope and check if the set is empty
           poly = mptPolytope(A_,b_);
           res = isempty(poly);
       end
              
   end 
end

%------------- END OF CODE --------------