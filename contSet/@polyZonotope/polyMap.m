function res = polyMap(pZ,coeff,E)
% polyMap - computes the polynomial map of a polyZonotope
%
% Syntax:
%    pZ = polyMap(pZ,coeff)
%    pZ = polyMap(pZ,coeff,E)
%
% Inputs:
%    pZ - polyZonotope object
%    coeff - coefficient matrix for the polynomial map 
%            p(x) = sum(coeff(:,i).*prod(x.^E,1),2)
%    E - exponent matrix for the polynomial map 
%             p(x) = sum(coeff(:,i).*prod(x.^E,1),2)
%
% Outputs:
%    pZ - polyZonotope object
%
% Example:
%    % polynomial zonotope
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    
%    % polynomial map [x(2)^3 + 2*x(1)*x(2)^2; -x(2)^3]
%    coeff = [1 2;-1 0];
%    E = [0 1;3 2];
%
%    % compute polynomial map
%    res = polyMap(pZ,coeff,E)
%
%    % visualization
%    plot(res);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, polyZonotope/cubMap

% Authors:       Niklas Kochdumper
% Written:       23-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check inpus
    inputArgsCheck({ ...
            {coeff, 'att', 'double', {'finite', 'matrix'}};
            {E, 'att', 'double', {'finite', 'matrix'}};
        })

    if size(E,1) ~= dim(pZ) || size(E,2) ~= size(coeff,2)
        throw(CORAerror('CORA:wrongValue', "second", ...
                 'dimension must match set and coefficient matrix.'));
    end

    % extract constant inputs
    n = dim(pZ);
    summands = {};
    c = zeros(size(coeff,1),1);
    tmp = sum(E,1);
    ind = find(tmp == 0);

    if ~isempty(ind)
        c = sum(coeff(:,ind));
        E(:,ind) = [];
        coeff(:,ind) = [];
    end

    % loop over all monomials
    for i = 1:size(E,2)
        
        % split monomial into a sequence of quadratic maps
        E_ = [];

        for j = 1:size(E,1)
            E_ = [E_, repelem(j,E(j,i))];
        end

        % compute results for all quadratic maps
        if mod(length(E_),2) == 1
            list = {project(pZ,E_(end))};
        else
            list = {};
        end

        for j = 2:2:length(E_)
            Q = zeros(n);
            Q(E_(j-1),E_(j)) = 1;
            list{end+1} = quadMap(pZ,{Q});
        end

        % compute product of all quadratic maps
        res = list{1};

        for j = 2:length(list)
            res = quadMap(res,list{j},{1},true);
        end

        % construct final set for this monomial
        summands{end+1} = coeff(:,i) * res;
    end

    % combine all monomials
    if ~isempty(summands)
        res = summands{1} + c;
    
        for i = 2:length(summands)
            res = exactPlus(res,summands{i});
        end
    else
        res = 0*pZ + c;
    end
end

% ------------------------------ END OF CODE ------------------------------
