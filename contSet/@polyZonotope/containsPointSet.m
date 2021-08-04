function res = containsPointSet(pZ,points,varargin)
% containsPointSet - checks if a point set is fully enclosed by a tight
%                    over-approximation of a polynomial zonotope
%
% Syntax:  
%    result = containsPointSet(pZ,points)
%    result = containsPointSet(pZ,points,splits,order)
%
% Inputs:
%    pZ - polyZonotope object (n dimensional)
%    points - matrix containing the point cloud (dimension: [n,M])
%    splits - number of splits for refinement (optional)
%    order - zonotope order for the linear zonotopes that over-approximate
%            the splitted polynomial zonotopes
%
% Outputs:
%    result - 1/0 if point is inside the zonotope or not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/containsPoint

% Author:       Niklas Kochdumper
% Written:      29-March-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values
splits = 4;
order = 10;

% parse input arguments
if nargin > 3 && ~isempty(varargin{1})
   splits = varargin{1};
end
if nargin == 4 
   order = varargin{2};
end

% split the reduced polynomial zonotope to obtain a better
% over-approximation with linear zonotopes
pZsplit{1} = pZ;

for i=1:splits
    qZnew = [];
    for j=1:length(pZsplit)
        try
            res = splitLongestGen(pZsplit{j});
        catch
            res = splitLongestIndepGen(pZsplit{j}); 
        end
        qZnew{end+1} = res{1};
        qZnew{end+1} = res{2};
    end
    pZsplit = qZnew;
end

% over-approximate the splitted polynomial zonotopes with linear zonotopes
for i = 1:length(pZsplit)
    temp = zonotope(pZsplit{i});
    temp = reduce(temp,'girard',order);
    qZsplit{i} = halfspace(temp);
end

temp = zonotope(pZ);
temp = reduce(temp,'girard',order);
qZsplit{end+1} = halfspace(temp);

% check if all points of the original polynomial zonotope are located
% indside the reduced polynomial zonotope
res = 1;

for i = 1:size(points,2)
   resTemp = 0;
   for j = 1:length(qZsplit)
      if in(qZsplit{j},points(:,i),1e-10)
         resTemp = 1;
         break;
      end
   end
   
   if resTemp == 0
      res = 0;
      break;
   end
end


% Auxiliary Functions -----------------------------------------------------

function splits = splitLongestIndepGen(pZ)
% split a polynomial zonotope at the longes independent generator
    
    % determine longest independent generator
    len = sum(pZ.Grest.^2,1);
    [~,ind] = max(len);
    ind = ind(1);
    g = pZ.Grest(:,ind);
    
    % split the polynomial zonotope
    pZ1 = pZ;
    pZ2 = pZ;
    
    pZ1.Grest(:,ind) = 0.5 * g;
    pZ1.c = pZ.c + 0.5 * g;
    
    pZ2.Grest(:,ind) = 0.5 * g;
    pZ2.c = pZ.c - 0.5 * g;

    splits{1} = pZ1;
    splits{2} = pZ2;

%------------- END OF CODE --------------