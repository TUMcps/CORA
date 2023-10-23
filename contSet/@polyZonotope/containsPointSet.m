function res = containsPointSet(pZ,points,varargin)
% containsPointSet - checks if a point set is fully enclosed by a tight
%    over-approximation of a polynomial zonotope
%
% Syntax:
%    res = containsPointSet(pZ,points)
%    res = containsPointSet(pZ,points,splits,order)
%
% Inputs:
%    pZ - polyZonotope object (n dimensional)
%    points - matrix containing the point cloud (dimension: [n,M])
%    splits - number of splits for refinement (optional)
%    order - zonotope order for the linear zonotopes that over-approximate
%            the splitted polynomial zonotopes
%
% Outputs:
%    res - true/false (over all points)
%
% Example:
%    pZ = polyZonotope([1;2],[1 -2 1; 2 3 1],[0;0],[1 0 2;0 1 1]);
%    points = [0 0 1; 1 2 4];
%    res = containsPointSet(pZ,points)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/containsPoint

% Authors:       Niklas Kochdumper
% Written:       29-March-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[splits,order] = setDefaultValues({4,10},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {points,'att','numeric'};
                {splits,'att','numeric','nonnan'};
                {order,'att','numeric','nonnan'}});

% split the reduced polynomial zonotope to obtain a better
% over-approximation with linear zonotopes
pZsplit{1} = pZ;

for i=1:splits
    qZnew = [];
    for j=1:length(pZsplit)
        try
            res = splitLongestGen(pZsplit{j});
        catch
            res = aux_splitLongestIndepGen(pZsplit{j}); 
        end
        qZnew{end+1} = res{1};
        qZnew{end+1} = res{2};
    end
    pZsplit = qZnew;
end

% over-approximate the splitted polynomial zonotopes with linear zonotopes
qZsplit = cell(length(pZsplit),1);
for i = 1:length(pZsplit)
    temp = zonotope(pZsplit{i});
    temp = reduce(temp,'girard',order);
    qZsplit{i} = halfspace(temp);
end

temp = zonotope(pZ);
temp = reduce(temp,'girard',order);
qZsplit{end+1} = halfspace(temp);

% check if all points of the original polynomial zonotope are located
% inside the reduced polynomial zonotope

% init container for point-in-splitZonotope
resTemp = false(length(qZsplit),size(points,2));
for j=1:length(qZsplit)
    resTemp(j,:) = contains_(qZsplit{j},points,'approx',1e-10);
end
% unify results: a point has to be in at least one of the split zonotopes
res = all(any(resTemp,1));

% for i = 1:size(points,2)
%    resTemp = false;
%    for j = 1:length(qZsplit)
%       if contains_(qZsplit{j},points(:,i),'approx',1e-10)
%          resTemp = true;
%          break;
%       end
%    end
%    
%    if ~resTemp
%       res = false;
%       break;
%    end
% end

end


% Auxiliary functions -----------------------------------------------------

function splits = aux_splitLongestIndepGen(pZ)
% split a polynomial zonotope at the longes independent generator
    
    % determine longest independent generator
    len = sum(pZ.GI.^2,1);
    [~,ind] = max(len);
    ind = ind(1);
    g = pZ.GI(:,ind);
    
    % split the polynomial zonotope
    pZ1 = pZ;
    pZ2 = pZ;
    
    pZ1.GI(:,ind) = 0.5 * g;
    pZ1.c = pZ.c + 0.5 * g;
    
    pZ2.GI(:,ind) = 0.5 * g;
    pZ2.c = pZ.c - 0.5 * g;

    splits{1} = pZ1;
    splits{2} = pZ2;

end

% ------------------------------ END OF CODE ------------------------------
