function val = supportFuncSplit(pZ,dir,type, varargin)
% supportFuncSplit - computes the support function of a (list of)
%    polyZonotopes in the given direction
%
% Syntax:
%    val = supportFuncSplit(pZ,dir,type)
%    val = supportFuncSplit(pZ,dir,type,splits)
%
% Inputs:
%    pZ - polyZonotope object or cell of polyZonotope objects
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    split - number of splits that are performed to calculate the bounds
%
% Outputs:
%    val - interval object specifying the upper and lower bound along the
%          direction
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/supportFunc

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      29-July-2018
% Last update:  17-October-2022  (NK: improve 'split' method)
%               06-December-2022 (TL: fix: 'split' considers splitted sets)
%               09-December-2022 (TL: speed up computation)
% Last revision: ---

%------------- BEGIN CODE --------------

splits = setDefaultValues({8}, varargin{:});

inputArgsCheck({ ...
    {pZ, 'att', {'cell', 'polyZonotope'}}; ...
    {dir,'att','numeric','vector'}; ...
    {type,'str',{'lower','upper','range'}}; ...
    {splits, 'att', 'numeric', {'integer', 'positive'}}
});

 % handle different types
if strcmp(type,'lower')
   val = -supportFuncSplit(pZ,-dir,'upper',splits);
   return;
elseif strcmp(type,'range')
   up = supportFuncSplit(pZ,dir,'upper',splits);
   low = -supportFuncSplit(pZ,-dir,'upper',splits);
   val = interval(low,up);
   return;
end

if isa(pZ, 'cell')
    pZsplit = pZ;
else
    pZsplit = {pZ};
end

% project the polynomial zonotope onto the direction
for i=1:length(pZsplit)
    pZsplit{i} = dir' * pZsplit{i};
end

% goal is to find a tight upper bound of all sets in pZsplit

% determine lower bound of upper bound
% is used to exclude entire sets with no need to split them further
minUpperBound = -inf;

% the result will be the smallest upper bound of the upper bound
maxUpperBound = -inf; % result for empty set

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape

for i = 0:splits
    % preinit to avoid copying
    qZnew = cell(2*length(pZsplit), 1);
    c = 0; % counter

    % reset to only consider leftover splitted subsets
    maxUpperBound = -inf;

    for j = 1:length(pZsplit) 
        if i == 0
            % check input without splitting (efficient for zonotopic input)
            res = pZsplit(j);
        else
            res = splitLongestGen(pZsplit{j});
        end
        
        for k = 1:length(res)
            res_k = res{k};
            
            % compute support function for enclosing zonotope
            [val_k,~,alpha] = supportFunc(zonotope(res_k),1);
            
            % update upper and lower bound
            maxUpperBound = max(maxUpperBound, val_k);
            
            if val_k >= minUpperBound
                % update min upper bound by 'most critical' point
                % aka largest point in zonotope subset

                % exract zonotopic generators from expMat
                ind1 = sum(res_k.expMat,1) == 1;
                ind2 = sum(res_k.expMat(:,ind1),2) == 1;
                alpha_ = zeros(size(res_k.expMat,1),1);
                alpha_(ind2) = alpha(ind1);
                
                % use result from zonotope supportFunc
                x_min = res_k.c + ...
                    sum(res_k.G .* prod(alpha_.^res_k.expMat,1)); 
                
                if ~isempty(res_k.Grest)
                    % same for Grest
                    beta = alpha(size(res_k.expMat,2)+1:end);
                    x_min = x_min + res_k.Grest*beta;
                end

                minUpperBound = max(minUpperBound, x_min);

                if withinTol(x_min, val_k)
                    % found exact upper bound for current set
                    continue;
                end
                
                % TODO other min upper bound? 
                % update min upper bound by largest random point (slow)
                % x_min = max(res_k.randPoint(11));
                % minUpperBound = max(minUpperBound, x_min);

                % add new set to queue
                c = c + 1;
                qZnew{c} = res_k; 
            end
        end
    end

    if withinTol(minUpperBound, maxUpperBound)
        % exact upper bound is found
        val = maxUpperBound;
        return
    end
    
    % update remaining splitted sets
    pZsplit = qZnew(1:c);
end

% return result
val = maxUpperBound;
end

%------------- END OF CODE --------------