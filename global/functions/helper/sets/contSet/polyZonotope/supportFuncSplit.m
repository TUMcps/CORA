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

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       29-July-2018
% Last update:   17-October-2022 (NK, improve 'split' method)
%                06-December-2022 (TL, fix: 'split' considers splitted sets)
%                09-December-2022 (TL, speed up computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

splits = setDefaultValues({8}, varargin);

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

% for optimizations, we also store the largest exact upper bound
maxExactUpperBound = -inf;

% split the polynomial zonotope multiple times to obtain a better 
% over-approximation of the real shape

for i = 0:splits
    % preinit to avoid copying
    qZnew = cell(2*length(pZsplit), 1);
    c = 0; % counter

    % reset to only consider leftover splitted subsets
    maxUpperBound = maxExactUpperBound;

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
            [max_k,~,alpha] = supportFunc(zonotope(res_k),1);
            
            % update upper and lower bound
            maxUpperBound = max(maxUpperBound, max_k);
            
            if max_k >= minUpperBound
                % update min upper bound by 'most critical' point
                % aka largest point in zonotope subset

                % exract zonotopic generators from E
                ind1 = sum(res_k.E,1) == 1;
                ind2 = sum(res_k.E(:,ind1),2) == 1;
                alpha_ = zeros(size(res_k.E,1),1);
                alpha_(ind2) = alpha(ind1);
                
                % use result from zonotope supportFunc
                minMax_k = res_k.c + ...
                    sum(res_k.G .* prod(alpha_.^res_k.E,1)); 
                
                if ~isempty(res_k.GI)
                    % same for GI
                    beta = alpha(size(res_k.E,2)+1:end);
                    minMax_k = minMax_k + res_k.GI*beta;
                end

                minUpperBound = max(minUpperBound, minMax_k);

                if withinTol(minMax_k, max_k)
                    % found exact upper bound for current set
                    maxExactUpperBound = max(maxExactUpperBound, max_k);
                    continue;
                end
                
                % TODO other min upper bound? 
                % update min upper bound by largest random point (slow)
                % minMax_k = max(res_k.randPoint(11));
                % minUpperBound = max(minUpperBound, minMax_k);

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

% ------------------------------ END OF CODE ------------------------------
