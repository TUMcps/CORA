function val = boundDir(obj,dir,type)
% boundDir - Calculate the upper or lower bound of a constrained zonotope
%            object along a certain direction
%
% Syntax:  
%    val = boundDir(obj,dir,type)
%
% Inputs:
%    obj - c-zonotope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % object properties  
    n = size(obj.A, 2);
    A = obj.A;
    b = obj.b;

    % ksi in [-1, 1]
    lb = -ones(n,1);
    ub = ones(n,1);
    
    % project the zonotope onto the direction
    f = dir' * obj.Z(:,2:end);
    
    % objective function equals zero => function value zero (linprog would
    % fail to find the solution in this case)
    if ~any(f)
        
        fval = 0;    
        
    else   
        
        % linear program options
        options = optimoptions('linprog','Algorithm','dual-simplex', 'display','off');

        % determine lower bound along the specified direction
        try
            if strcmp(type,'lower')
                [~, fval,exitflag] = linprog(f',[],[],A,b,lb,ub,options);
            else
                [~, fval,exitflag] = linprog(-f',[],[],A,b,lb,ub,options);
                fval = -fval;
            end
        catch
           exitflag = 0;
           fval = [];
        end

        % check if a solution could be found
        if exitflag <= 0
            
            % check if the objective function is identical to the normal vector
            % of a constraint -> only one solution
            for i = 1:size(A,1)
               if sum(abs(f-A(i,:))) < 1e-15 
                  fval = b(i); 
                  break;
               elseif sum(abs(f+A(i,:))) < 1e-15
                  fval = -b(i);
                  break;
               end
            end
            
            % try out different algorithms if "dual-simplex" failed
            if isempty(fval)
               
                % linear program options
                options = optimoptions('linprog','Algorithm','interior-point', 'display','off');

                % determine lower bound along the specified direction
                if strcmp(type,'lower')
                    [~, fval,exitflag] = linprog(f',[],[],A,b,lb,ub,options);
                else
                    [~, fval,exitflag] = linprog(-f',[],[],A,b,lb,ub,options);
                    fval = -fval;
                end  
                
                % error message if still no solution
                if exitflag <= 0
                   error('conZonotope/boundDir: linear program failed to find a solution!'); 
                end
            end
        end
    end
    
    % calculate bound by adding the zonotope center
    val = dir' * obj.Z(:,1) + fval;

end

%------------- END OF CODE --------------