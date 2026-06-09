function [r,grad] = gradientRobustness(phi,varargin)
% gradientRobustness - compute the gradient of a temporal logic formula 
%    with respect to the states of a trajectory (see [1])
%
% Syntax:
%    [r,grad] = gradientRobustness(phi,traj)
%    [r,grad] = gradientRobustness(phi,x,t)
%
% Inputs:
%    phi - temporal logic formula (class stl)
%    traj - simulated trajectory (class trajectory)
%    x - states of the trace (dimensions: [n,m])
%    t - times of the trace (dimensions: [1,m])
%
% Outputs:
%    r - robustness of the STL formula for the current trajectory
%    grad - gradient of the robustness for the STL formula with repect to
%           the states x of the trajectory
%
% Example: 
%    x = stl('x',2);
%    phi = finally(globally(x(1) > 15,interval(0,0.2)),interval(0,3));
%
%    sys = linearSysDT([0.72 0.36; -0.18 1.08],0.2);
%
%    simOpts.x0 = [-10;10];
%    simOpts.tFinal = 5;
%    [t,x] = simulate(sys,simOpts);
%
%    [r,grad] = gradientRobustness(phi,x,t)
%
% References: 
%   [1] Y. V. Pant and et al, "Smooth Operator: Control using the Smooth 
%       Robustness of Temporal Logic", CCTA 2017   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       11-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    if nargin == 2
        traj = varargin{1};
        x = traj(1).x; t = traj(1).t;
    else
        x = varargin{1}; t = varargin{2};
    end

    % bring temporal logic formula to correct format and extract predicates
    [phi,~,sets] = aux_preprocessTemporalLogic(phi);

    % precompute gradient for all predicates
    r = zeros(length(sets),size(x,2));
    grad = cell(length(sets),1);

    for i = 1:length(sets)

        grad{i} = zeros(size(x'));

        for j = 1:size(x,2)
            [r(i,j),grad{i}(j,:)] = aux_gradientPredicate(x(:,j),sets{i});
        end
    end

    % compute robustness and gradient for the STL formula
    [r,grad] = aux_gradientTemporalLogic(phi,r,grad,t);

    r = r(1); grad = grad{1}';
end


% Auxiliary functions -----------------------------------------------------

function [r,grad] = aux_gradientTemporalLogic(phi,r,grad,time)
% recursive function to compute the gradient of the robustness of a 
% temporal logic formula

    if ~phi.temporal

        r = r(phi.id,:);
        gradTmp = grad{phi.id};

        grad = cell(1,length(r));

        for i = 1:length(r)
            grad{i} = zeros(length(r),size(gradTmp,2));
            grad{i}(i,:) = gradTmp(i,:);
        end

    elseif strcmp(phi.type,'&') % ---
        % compute robustness of each hs
        [r1,grad1] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);
        [r2,grad2] = aux_gradientTemporalLogic(phi.rhs,r,grad,time);

        r = zeros(size(r1)); grad = cell(1,length(r1));

        for i = 1:length(r1)
            [r(i),gradTmp] = aux_gradientSoftmin([r1(i);r2(i)]);
            grad{i} = gradTmp(1)*grad1{i} + gradTmp(2)*grad2{i};
        end

    elseif strcmp(phi.type,'|') % ---
        % compute robustness of each hs
        [r1,grad1] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);
        [r2,grad2] = aux_gradientTemporalLogic(phi.rhs,r,grad,time);

        r = zeros(size(r1)); grad = cell(1,length(r1));

        for i = 1:length(r1)
            [r(i),gradTmp] = aux_gradientSoftmax([r1(i);r2(i)]);
            grad{i} = gradTmp(1)*grad1{i} + gradTmp(2)*grad2{i};
        end

    elseif strcmp(phi.type,'next') % ---

        [r,grad] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);

        index = find(time >= phi.from);
        index = intersect(index,1:length(r));

        len = length(r)-length(index);
        r = [r(index), -inf*ones(1,len)];
        grad = [grad(index),repmat({zeros(size(grad{1}))},[1,len])];

    elseif strcmp(phi.type,'finally') % --- 

        [r_,grad_] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r_));
        grad = repmat({zeros(size(grad_{1}))},[1,length(r)]);

        while ~isempty(index) && index(1) <= length(r)
            index = index(index <= length(r));
            [r(cnt),gradTmp] = aux_gradientSoftmax(r_(index)');
            for i = 1:length(gradTmp)
                grad{cnt} = grad{cnt} + grad_{index(i)} * gradTmp(i);
            end
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'globally') % ---

        [r_,grad_] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r_));
        grad = repmat({zeros(size(grad_{1}))},[1,length(r)]);

        while ~isempty(index) && index(1) <= length(r)
            index = index(index <= length(r));
            [r(cnt),gradTmp] = aux_gradientSoftmin(r_(index)');
            for i = 1:length(gradTmp)
                grad{cnt} = grad{cnt} + grad_{index(i)} * gradTmp(i);
            end
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'until') % ---

        % compute robustness of each hs
        [r1,grad1] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);
        [r2,grad2] = aux_gradientTemporalLogic(phi.rhs,r,grad,time);
        
        % find indices
        index = find(time >= phi.from & time <= phi.to);

        % init
        cnt = 1; 
        r = -inf * ones(size(r1));
        grad = repmat({zeros(size(grad1{1}))},[1,length(r)]);
    
        % check each index
        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));
            r_ = -inf * ones(1,length(index));
            g_ = cell(1,length(index));

            % get max-min robustness
            for i = 1:length(index)
                [r_(i),gradTmp] = aux_gradientSoftmin([r2(index(i))';r1(cnt:index(i))']);
                g_{i} = grad2{index(i)}*gradTmp(1);
                for j = 2:length(gradTmp)
                    g_{i} = g_{i} + grad1{cnt+(j-2)}*gradTmp(j);
                end
            end

            [r(cnt),gradTmp] = aux_gradientSoftmax(r_');

            for i = 1:length(index)
                grad{cnt} = grad{cnt} + gradTmp(i)*g_{i};
            end
    
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'release') % ---

        % compute robustness of each hs
        [r1,grad1] = aux_gradientTemporalLogic(phi.lhs,r,grad,time);
        [r2,grad2] = aux_gradientTemporalLogic(phi.rhs,r,grad,time);
        
        % find indices
        index = find(time >= phi.from & time <= phi.to);

        % init 
        cnt = 1; 
        r = -inf * ones(size(r1));
        grad = repmat({zeros(size(grad1{1}))},[1,length(r)]);
        r1 = -r1; r2 = -r2;
    
        % check each index
        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));

            r_ = -inf * ones(1,length(index));
            g_ = cell(1,length(index));
    
            for i = 1:length(index)
                % get max-min robustness
                [r_(i),gradTmp] = aux_gradientSoftmin([r2(index(i))';r1(cnt:index(i))']);
                g_{i} = -grad2{index(i)}*gradTmp(1);
                for j = 2:length(gradTmp)
                    g_{i} = g_{i} - grad1{cnt+(j-2)}*gradTmp(j);
                end
            end

            [r(cnt),gradTmp] = aux_gradientSoftmax(r_');

            for i = 1:length(index)
                grad{cnt} = grad{cnt} + gradTmp(i)*g_{i};
            end
    
            cnt = cnt + 1; index = index + 1;
        end

        % negate result?
        r(~isinf(r)) = -r(~isinf(r));
    end
end

function [phi,pred,sets] = aux_preprocessTemporalLogic(phi)
% preprocess temporal logic formula

    % convert to negation normal form
    phi = negationNormalForm(phi);

    % assign unique identifiers to all predicates
    [phi,pred] = assignIdentifiers(phi);

    % convert the regions defined by the predicates to sets
    sets = cell(size(pred));

    for i = 1:length(pred)

        % convert to a union of safe sets
        eq = disjunctiveNormalForm(pred{i});
        clauses = getClauses(eq,'dnf');

        if length(clauses) == 1                 % single safe set

            safeSet = convert2set(clauses{1});
            sets{i} = aux_reverseInequalityConstraints(safeSet);

        else                                    % union of safe sets

            list = cell(length(clauses),1);

            for j = 1:length(clauses)
                list{j} = convert2set(clauses{j});
            end

            % convert to a union of unsafe sets
            sets{i} = aux_safe2unsafe(list);
        end
    end
end

function list = aux_safe2unsafe(sets)
% convert a safe set defined by the union of multiple sets to an
% equivalent union of unsafe sets

    % reverse first constraint
    list = aux_reverseInequalityConstraints(sets{1});

    for i = 2:length(sets)

        % reverse next constraint
        nextConstReverse = aux_reverseInequalityConstraints(sets{i});

        % go through all combinations
        list_ = {};
        for j = 1:length(nextConstReverse)
            for k = 1:length(list)
                if isa(list{k},'levelSet') || isa(nextConstReverse{j},'levelSet') || ...
                        isIntersecting_(list{k},nextConstReverse{j},'exact',1e-8)
                    % compute intersection
                    list_{end+1} = and_(list{k},nextConstReverse{j},'exact');
                end
            end
        end

        % update list
        list = list_;
    end
end

function res = aux_reverseInequalityConstraints(S)
% get a list of reversed inequality constraints for a given set

    res = {};

    if isa(S,'levelSet')
        compOp = S.compOp;

        if ~iscell(compOp)
           compOp = {compOp};
        end

        for i = 1:size(S.eq,1)
            res{end+1} = levelSet(-S.eq(i),S.vars,compOp{i});
        end

    else
        % convert to polytope
        poly = polytope(S);
        for i = 1:length(poly.b)
            res{end+1} = ~polytope(poly.A(i,:),poly.b(i));
            res{end} = normalizeConstraints(res{end},'A');
        end
    end
end

function [val,grad] = aux_gradientSoftmax(x)
% compute the gradient of the softmax function 

    % compute softmax
    expX = exp(x);
    s = expX./sum(expX);
    val = s'*x;

    % compute Jacobian function according to https://eli.thegreenplace.net/
    % 2016/the-softmax-function-and-its-derivative/
    J = -s*s';
    J = J - diag(diag(J)) + diag(s.*(1-s));

    % compute gradient for sum_i s(i)*x(i)
    grad = sum(J.*x + diag(s),1);
end

function [val,grad] = aux_gradientSoftmin(x)
% compute the gradient of the softmin function 

    % compute softmax
    expNegX = exp(-x);
    s = expNegX./sum(expNegX);
    val = s'*x;

    % compute Jacobian function according to https://eli.thegreenplace.net/
    % 2016/the-softmax-function-and-its-derivative/
    J = -s*s';
    J = J - diag(diag(J)) + diag(s.*(1-s));

    % compute gradient for sum_i s(i)*x(i)
    grad = sum(-J.*x + diag(s),1);
end

function [val,grad] = aux_gradientPredicate(x,sets)
% compute the gradient for a predicate (= distance to an unsafe set)

    % single unsafe set
    if length(sets) == 1 
        
        % only one halfspace constraint
        if size(sets{1}.A) == 1

            val = sets{1}.A*x - sets{1}.b;
            grad = sets{1}.A;

        % multiple halfspace constraints
        else
            halfspaceVals = sets{1}.A*x - sets{1}.b;
            [val,grad] = aux_gradientSoftmax(halfspaceVals);
            grad = sum(diag(grad)*sets{1}.A,1);
        end

    % multiple unsafe sets
    else

        % loop over all unsafe sets
        gradTmp = zeros(length(sets),length(x)); predVals = zeros(length(sets));

        for i = 1:lenght(sets)
            [predVals(i),gradTmp(i,:)] = aux_gradientPredicate(x,sets(i));
        end

        % take the minimum distance to any unsafe set
        [val,grad] = aux_gradientSoftmin(predVals);
        grad = sum(diag(grad)*gradTmp,1);
    end
end

% ------------------------------ END OF CODE ------------------------------
