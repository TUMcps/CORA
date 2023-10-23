function list = generateRandomPredicates(n,nrPred,dom)
% generateRandomPredicates - generate random predicates for an STL formula
%
% Syntax:
%    list = generateRandomPredicates(n,nrOps,dom)
%
% Inputs:
%    n - dimension of the state space
%    nrPred - number of predicates
%    dom - domain where predicates flip between true/false (class interval)
%
% Outputs:
%    list - list of random predicates (cell-array of class stl)
%
% Example: 
%    list = generateRandomPredicates(2,3,interval([-1;-1],[1;1]));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       15-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % generate predicates
    x = stl('x',n);
    list = cell(nrPred,1);
    operators = {'<','<=','>','>='};

    for i = 1:nrPred

        % choose random operator
        op = operators{randi(4)};

        % halfspace constraint (single variable or all variables)
        if rand() > 0.5

            % choose variable
            ind = randi(n);

            % generate offset
            if representsa_(dom,'emptySet',eps)
                offset = -5 + 10*rand();
            else
                offset = randPoint(dom(ind));
            end

            % generate predicate
            eval(['list{i} = x(ind) ',op,num2str(offset),';']);

        else

            % choose halfspace normal vector
            c = -2 + 4*rand(n,1);

            % choose offset
            if representsa_(dom,'emptySet',eps)
                offset = -2 + 4*rand();
            else
                offset = randPoint(c'*zonotope(dom));
            end

            % generate predicate
            rhs = c(1) * x(1);

            for j = 2:n
                rhs = rhs + c(j) * x(j);
            end
            
            eval(['list{i} = rhs',op,num2str(offset),';']);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
