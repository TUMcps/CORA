function I = tightenDomain(ls,S,varargin)
% tightenDomain - Computes a tight interval enclosure of the intersection 
%    between the level set and the set S as described in Sec. 4.2 in [1];
%    this operation only works for level sets with '=='
%
% Syntax:
%    I = tightenDomain(ls,S)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object
%
% Outputs:
%    I - interval over-approximation of the intersection
%
% Example: 
%    syms x y
%    eq = y*sin(x)^2 + x*y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    Z = zonotope([0 0.5 0 0.5;2 0 0.5 0.5]);
% 
%    I1 = tightenDomain(ls,Z,'forwardBackward');
%    I2 = tightenDomain(ls,Z,'split');
% 
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(Z,[1,2],'FaceColor','g');
%    plot(ls,[1,2],'b');
%    plot(I1,[1,2],'r');
%    plot(I2,[1,2],'k');
%
% References:
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and

% Authors:       Niklas Kochdumper
% Written:       25-July-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% operations only if '=='
if ~strcmp(ls.compOp,'==')
    throw(CORAerror('CORA:noops',ls,S));
end

% default values
maxIter = 8;
type = setDefaultValues({'forwardBackward'},varargin);

% check input arguments
inputArgsCheck({{ls,'att','levelSet'};
                {S,'att','contSet'};
                {type,'str',{'forwardBackward','split','taylor','linear'}}});

% enclose set by interval
I = interval(S);

% choose contractor algorithm 
if strcmp(type,'split')
    I = aux_contractorSplit(ls,I,maxIter);
elseif strcmp(type,'taylor')
    I = aux_contractorTaylor(ls,I,maxIter);
elseif strcmp(type,'linear')
    I = aux_contractorLinear(ls,I,maxIter);
elseif strcmp(type,'forwardBackward')
    I = aux_contractorForwardBackward(ls,I);
end
    
end
   

% Auxiliary functions -----------------------------------------------------

function int = aux_contractorForwardBackward(obj,int)

    % compute syntax tree nodes for all variables
    vars = [];
    
    for i = 1:length(int)
        vars = [vars;syntaxTree(int(i),i)];
    end

    % forward iteration: compute syntax tree
    synTree = obj.funHan(vars);
    
    % backward iteration
    int = backpropagation(synTree,interval(0,0),int);

end

function int = aux_contractorLinear(obj,int,maxIter)

    % loop over all iterations
    for i = 1:maxIter
       
        % evaluate gradient 
        m = center(int);
        k = obj.der.grad(int);
        f = obj.funHan(m);
        
        sup_ = supremum(k);
        inf_ = infimum(k);
        
        % loop over all variables
        for j = 1:length(k)
           
            if sign(sup_(j)) == sign(inf_(j))
                k_ = k;
                k_(j) = interval(0,0);
                temp = (-f + transpose(k)*m - transpose(k_)*int)/k(j);
                int(j) = and_(int(j),temp,'exact');
            end
        end
    end
end

function int = aux_contractorTaylor(obj,int,maxIter)

    infi = infimum(int);
    sup = supremum(int);

    % loop over all iterations
    for i = 1:maxIter
            
        % loop over all variables
        for j = 1:length(obj.solved)

            % check if level set equation contains the current variable
            if obj.solved{j}.contained

                % check if equaiton is solvable for the current variable
                if obj.solved{j}.solvable

                    [infi,sup,fail] = aux_scaleVarSolvable(obj,infi,sup,j);

                end

                % not solvable or evaluation failed
                if ~obj.solved{j}.solvable || fail
                    [infi,sup] = aux_scaleVarUnsolvable(obj,infi,sup,j);
                end
            end
        end            
    end
    
    % assign final interval
    int = interval(infi,sup);
end
    
function int = aux_contractorSplit(obj,int,maxIter)   
    
    intList{1} = int;
    
    % loop over all iterations
    for i = 1:maxIter
        
        intNew = {};
        
        % loop over all parallel intervals
        for h = 1:length(intList)
            
            infi = infimum(intList{h});
            sup = supremum(intList{h});
            
            % loop over all variables
            for j = 1:length(obj.solved)

                % check if level set equation contains the current variable
                if obj.solved{j}.contained

                    % check if equaiton is solvable for the current variable
                    if obj.solved{j}.solvable

                        [infi,sup,fail] = aux_scaleVarSolvable(obj,infi,sup,j);
                        
                    end
                        
                    % not solvable or evaluation failed
                    if ~obj.solved{j}.solvable || fail

                        [infi,sup] = aux_scaleVarUnsolvable(obj,infi,sup,j);
                        
                        if infi(j) > sup(j)
                           break; 
                        end
                    end
                end
            end
            
            % split intervals
            if ~any(infi > sup)
                
                int_ = interval(infi,sup);
                vol_ = volume_(int_);

                if vol_ >= volume(intList{h})

                    % split along largest dimension
                    [~,ind] = max(rad(int_));
                    temp = split(int_,ind);
                    intNew{end+1} = temp{1};
                    intNew{end+1} = temp{2};
                   
                else
                   
                    intNew{end+1} = int_;       
                    
                end
            end
        end
        
        intList = intNew;
    end
    
    % compute union of all intervals
    if length(intList) >= 1
        int = intList{1};

        for i = 2:length(intList)
           int = int | intList{i}; 
        end
    else
       int = []; 
    end
end


% Auxiliary Functions -----------------------------------------------------

function [infi,sup,fail] = aux_scaleVarSolvable(obj,infi,sup,j)

    fail = 0;

    % loop over all possible solutions
    I = {};

    for k = 1:length(obj.solved{j}.funHan)

        funHan = obj.solved{j}.funHan{k}.eq;

        try
            intTemp = funHan(interval(infi,sup));

            if ~(infimum(intTemp) > sup(j) || ...
                 supremum(intTemp) < infi(j))

                I{end+1} = intTemp;
            end
        end
    end

    % compute union of all possible solutions
    if length(I) >= 1

        intTemp = I{1};

        for k = 2:length(I)
           intTemp = intTemp | I{k}; 
        end

        infi(j) = max(infimum(intTemp),infi(j));
        sup(j) = min(supremum(intTemp),sup(j));
    else
        fail = 1;
    end
end

function [infi,sup] = aux_scaleVarUnsolvable(obj,infi,sup,j)

    % compute Taylor series at linearization point
    int = interval(infi,sup);
    c = center(int);
    int_ = int - c;

    eq = obj.funHan(c);
    grad = obj.der.grad(c);
    hess = obj.der.hess(c);

    third = cell(length(c),1);
    for k = 1:length(third)
       han = obj.der.third{k};
       third{k} = han(int);
    end

    if grad(j) ~= 0
        
        % over-approximate higher order terms
        err = aux_quadEval(hess,int_);

        for k = 1:length(third)
            err = err + 1/6 * int_(k) * aux_quadEval(third{k},int_);
        end

        % update bounds
        ind = setdiff(1:length(c),j);
        intNew = (grad(j)*c(j) - eq - err - grad(ind)' * int_(ind))/grad(j);

        infi(j) = max(infimum(intNew),infi(j));
        sup(j) = min(supremum(intNew),sup(j));
    end

end

function res = aux_quadEval(Q,int)
% tight evaluation of a quadratic term using interval arithmetic

    res = interval(0,0);
        
    for k = 1:length(int)
        temp = interval(0,0);
        for l = k+1:length(int)
            temp = temp + 0.5*(Q(k,l) + Q(l,k)) * int(l);
        end
        res = res + 0.5 * Q(k,k) * int(k)^2 + temp * int(k);
    end
end

% ------------------------------ END OF CODE ------------------------------
