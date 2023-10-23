function res = contains_(ls,S,type,tol,varargin)
% contains_ - determines if a level set contains a set or a point
%
% Syntax:
%    res = contains_(ls,S)
%    res = contains_(ls,S,type)
%    res = contains_(ls,S,type,tol)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object or single point
%    type - 'exact' or 'approx'
%    tol - numerical tolerance (only for point in set containment)
%
% Outputs:
%    res - true/false
%
% Example: 
%    syms x y
%    eq = sin(x) + y;
%    ls = levelSet(eq,[x;y],'<=');
%
%    I1 = interval([0.7;-0.3],[1.3;0.3]);
%    I2 = interval([-1.3;-0.3],[-0.7;0.3]);
%
%    contains(ls,I1)
%    contains(ls,I2)
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I1,[1,2],'FaceColor','r');
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I2,[1,2],'FaceColor','g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, halfspace/contains_, conHyperplane/contains_

% Authors:       Niklas Kochdumper
% Written:       19-July-2019
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% set or single point 
if ~isnumeric(S)                                     % set

    % check type of level set
    if strcmp(ls.compOp,'==')
        throw(CORAerror('CORA:noops',ls,S));
    end

    % interval over-approximation
    I = interval(S); 

    % evaluate non-linear function with interval arithmetic
    I = ls.funHan(I);

    ub = supremum(I);

    % multiple inequality constraints or not
    if ~iscell(ls.compOp)

        % switch case for the different types of level sets
        if strcmp(ls.compOp,'<=')
            res = ub < 0 | withinTol(ub,0);
        else
            res = ub < 0;
        end

    else

        resVec = false(length(ls.compOp),1);

        % loop over all inequality constraints
        for i = 1:length(ls.compOp)

            if strcmp(ls.compOp{i},'<=')
                resVec(i) = ub(i) < 0 | withinTol(ub(i),0);
            else
                resVec(i) = ub(i) < 0;
            end
        end

        res = all(resVec);
    end

else                                                    % array of points

    % init resulting logical array
    res = false(1,size(S,2));

    for j=1:size(S,2)
    
        % evaluate nonlinear function
        val = ls.funHan(S);
        
         % multiple inequality constraints or not
        if ~iscell(ls.compOp)
    
            % switch case for the different types of level sets
            if strcmp(ls.compOp,'==')
                tmp = abs(val);
                res(j) = tmp < tol | withinTol(tmp,tol);
            elseif strcmp(ls.compOp,'<=')
                res(j) = val < tol | withinTol(val,tol);
            else
                res(j) = val < tol;
            end
    
        else
    
            resVec = false(length(ls.compOp),1);
    
            % loop over all inequality constraints
            for i = 1:length(ls.compOp)
    
                if strcmp(ls.compOp{i},'==')
                    tmp = abs(val(i));
                    resVec(i) = tmp < tol | withinTol(tmp,tol);
                elseif strcmp(ls.compOp{i},'<=')
                    resVec(i) = val(i) < tol | withinTol(val(i),tol);
                else
                    resVec(i) = val(i) < tol;
                end
            end
    
            res(j) = all(resVec);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
