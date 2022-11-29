function [val,x,ksi] = supportFunc(cZ,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a constrained 
%    zonotope along a certain direction
%
% Syntax:  
%    [val,x,ksi] = supportFunc(cZ,dir)
%    [val,x,ksi] = supportFunc(cZ,dir,type)
%
% Inputs:
%    cZ - conZonotope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constrained zonotope in the specified direction
%    x - support vector
%    ksi - factor values that correspond to the bound
%
% Example:
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    dir = [1;1];
%    val = supportFunc(cZ,dir);
%
%    figure; hold on; xlim([-3,1]); ylim([-3,4]);
%    plot(cZ);
%    plot(conHyperplane(dir,val),[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% pre-processing
[res,vars] = pre_supportFunc('conZonotope',cZ,dir,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = vars{1}; return
else
    cZ = vars{1}; dir = vars{2}; type = vars{3};
end

% range not supported
if strcmp(type,'range')
    throw(CORAerror('CORA:notSupported',type));
end


% constrained vs. non-constrained case
if isempty(cZ.A) || all(all(cZ.A == 0)) 
    
    % project zonotope onto the direction
    temp = dir'*cZ;
    inter = interval(temp);
    
    % determine upper or lower bound
    if strcmp(type,'lower')
       val = infimum(inter); 
       ksi = -sign(temp.Z(2:end))';
    else
       val = supremum(inter);
       ksi = sign(temp.Z(2:end))';
    end
    
else
    
    % object properties  
    n = size(cZ.A, 2);
    A = cZ.A;
    b = cZ.b;

    % ksi in [-1, 1]
    lb = -ones(n,1);
    ub = ones(n,1);

    % project the zonotope onto the direction
    f = dir' * cZ.Z(:,2:end);

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
                [ksi, fval,exitflag] = linprog(f',[],[],A,b,lb,ub,options);
            else
                [ksi, fval,exitflag] = linprog(-f',[],[],A,b,lb,ub,options);
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
                options = optimoptions('linprog','Algorithm','interior-point', ...
                                       'MaxIterations',10000,'display','off');

                % determine lower bound along the specified direction
                if strcmp(type,'lower')
                    [ksi,fval,exitflag] = linprog(f',[],[],A,b,lb,ub,options);
                else
                    [ksi,fval,exitflag] = linprog(-f',[],[],A,b,lb,ub,options);
                    fval = -fval;
                end  

                % error message if still no solution
                if exitflag <= 0
                    throw(CORAerror('CORA:solverIssue'));
                end
            end
        end
    end

    % calculate bound by adding the zonotope center
    val = dir' * cZ.Z(:,1) + fval;
end

% calculate support vector
if nargout >= 2
    x = cZ.Z(:,1) + cZ.Z(:,2:end)*ksi;
end

%------------- END OF CODE --------------