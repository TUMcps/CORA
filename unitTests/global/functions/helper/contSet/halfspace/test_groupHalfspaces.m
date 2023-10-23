function res = test_groupHalfspaces()
% test_groupHalfspaces - unit test function for groupHalfspaces
%
% Syntax:
%    res = test_groupHalfspaces()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: groupHalfspaces

% Authors:       Niklas Kochdumper
% Written:       24-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = true;

     % loop over different dimensions
     for j = 2:5
    
        % generate list of random halfspaces that are similar
        hs = halfspace.generateRandom('Dimension',j);
        
        c = interval(hs.c - 0.1*abs(hs.c), hs.c + 0.1*abs(hs.c));
        d = interval(hs.d - 0.1*abs(hs.d), hs.d + 0.1*abs(hs.d));
        
        list = cell(1 + randi(10),1);
        list{1} = hs;
        
        for i = 2:length(list)
            list{i} = halfspace(randPoint(c),randPoint(d));
        end
        
        % generate random domain that intersects the halfspaces
        dom = interval.generateRandom('Dimension',j);
        
        if ~isIntersecting(conHyperplane(hs),dom)
            n = norm(hs.c);
            tmp = interval((hs.c' / n) * zonotope(dom));
            dom = dom + (-center(dom)) + (hs.c / n) * randPoint(tmp); 
        end
        
        % group the halfpaces together
        hs_i = groupHalfspaces(list,dom,'inner');
        hs_o = groupHalfspaces(list,dom,'outer');
        
        % check the result for the inner-approximation
        tmp = hs_i & polytope(dom);
        
        if ~representsa(tmp,'emptySet')
            points = randPoint(tmp,10);
        
            for i = 1:length(list)
                if ~contains(list{i},points)
                    error('test_groupHalfspaces failed!');
                end
            end
        end
        
        % check the result for the outer-approximation
        points = [];
        
        for i = 1:length(list)
            tmp = list{i} & polytope(dom);
        
            if ~representsa(tmp,'emptySet')
                points = [points, randPoint(tmp,10)];
            end
        end
        
        if ~isempty(points)
            if ~contains(hs_o,points)
                error('test_groupHalfspaces failed!');
            end
        end
     end

% ------------------------------ END OF CODE ------------------------------
