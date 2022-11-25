function res = in(obj1,obj2,varargin)
% in - determines if obj2 is fully inside a levelSet obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%    res = in(obj1,p,tol)
%
% Inputs:
%    obj1 - levelSet object
%    obj2 - contSet object
%    p - single point
%    tol - numerical tolerance for point in set containment
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    syms x y
%    eq = sin(x) + y;
%    ls = levelSet(eq,[x;y],'<=');
%
%    int1 = interval([0.7;-0.3],[1.3;0.3]);
%    int2 = interval([-1.3;-0.3],[-0.7;0.3]);
%
%    in(ls,int1)
%    in(ls,int2)
%
%    figure
%    hold on
%    xlim([-1.5,1.5]);
%    ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(int1,[1,2],'r','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    xlim([-1.5,1.5]);
%    ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(int2,[1,2],'g','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/in, conHyperplane/in

% Author:       Niklas Kochdumper
% Written:      19-July-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % set or single point 
    if ~isnumeric(obj2)                                     % set

        % check type of level set
        if strcmp(obj1.compOp,'==')
           error('Operation not applicable for equality level sets!'); 
        end

        % interval over-approximation
        int = interval(obj2); 

        % evaluate non-linear function with interval arithmetic
        intRes = obj1.funHan(int);

        sup = supremum(intRes);

        % multiple inequality constraints or not
        if ~iscell(obj1.compOp)

            % switch case for the different types of level sets
            if strcmp(obj1.compOp,'<=')
                res = (sup <= 0);
            else
                res = (sup < 0);
            end

        else

            resVec = zeros(length(obj1.compOp),1);

            % loop over all inequality constraints
            for i = 1:length(obj1.compOp)

                if strcmp(obj1.compOp{i},'<=')
                    resVec(i) = (sup(i) <= 0);
                else
                    resVec(i) = (sup(i) < 0);
                end
            end

            res = all(resVec == 1);
        end
    
    else                                                    % single point
        
        % parse input arguments
        tol = 0;
        
        if nargin >= 3
            tol = varargin{1};
        end
        
        % evaluate nonlinear function
        val = obj1.funHan(obj2);
        
         % multiple inequality constraints or not
        if ~iscell(obj1.compOp)

            % switch case for the different types of level sets
            if strcmp(obj1.compOp,'==')
                res = (abs(val) <= tol);
            elseif strcmp(obj1.compOp,'<=')
                res = (val <= tol);
            else
                res = (val < tol);
            end

        else

            resVec = zeros(length(obj1.compOp),1);

            % loop over all inequality constraints
            for i = 1:length(obj1.compOp)

                if strcmp(obj1.compOp{i},'==')
                    resVec(i) = (abs(val(i)) <= tol);
                elseif strcmp(obj1.compOp{i},'<=')
                    resVec(i) = (val(i) <= tol);
                else
                    resVec(i) = (val(i) < tol);
                end
            end

            res = all(resVec == 1);
        end
    end

%------------- END OF CODE --------------