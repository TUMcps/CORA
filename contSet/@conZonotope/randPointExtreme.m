function p = randPointExtreme(cZ)
% randPointExtreme - generates a random extreme point inside a constrained 
%                    zonotope
%
% Syntax:  
%    p = randPointExtreme(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    p - random extreme point in R^n
%
% Example: 
%    cZ = conZonotope.generateRandom(2);
%    p = randPointExtreme(cZ);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(cZ.A)   

        % no constraints -> call superclass method
        p = randPointExtreme(zonotope(cZ.Z));

    else                

        % center constrained zonotope at origin
        c = center(cZ);
        cZ = cZ + (-c);
       
        % select random direction
        n = length(c);
        d = rand(n,1) - 0.5*ones(n,1);
        d = d./norm(d);
        
        % compute farest point in this direction that is still in set
        [~,x] = supportFunc(cZ,d);
        p = x + c;
    end
end

%------------- END OF CODE ----------