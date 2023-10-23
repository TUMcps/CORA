function val = hausdorffDist(pZ,points,varargin)
% hausdorffDist - Calculates an approximation of the Hausdorff distance
%    between a polynomial zonotope and a point cloud
%
% Syntax:
%    val = hausdorffDist(pZ,points)
%    val = hausdorffDist(pZ,points,splits)
%
% Inputs:
%    pZ - polyZonotope object (subset R^n)
%    points - points cloud provided as a matrix (dimension: [n,M])
%    splits - number of set splits that are performed to calculate a better
%             approximaion
%
% Outputs:
%    val - approximated Hausdorff distance
%
% Examples:
%    pZ = polyZonotope([0;0],[2 0 2;0 2 2],[0;0],[1 0 3;0 1 1]);
%    points = randPoint(pZ,10000);
%    dist = hausdorffDist(pZ,points,8)
% 
%    figure; hold on;
%    plot(points(1,:),points(2,:),'.k');
%    plot(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadZonotope/hausdorffDist

% Authors:       Niklas Kochdumper
% Written:       05-September-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    splits = setDefaultValues({4},varargin);

    % check input arguments
    inputArgsCheck({{pZ,'att','polyZonotope'};
                    {points,'att','numeric','matrix'};
                    {splits,'att','numeric','nonempty'}});
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ;

    for i=1:splits
        qZnew = [];
        for j=1:length(pZsplit)
            res = splitLongestGen(pZsplit{j});

            qZnew{end+1} = res{1};
            qZnew{end+1} = res{2};
        end
        pZsplit = qZnew;
    end
    
    % calculate the vertices of the zonotope over-approximations of the
    % splitted polynomial zonotopes
    vert = [];
    
    for i = 1:length(pZsplit)
        zono = zonotope(pZsplit{i});
        temp = vertices(zono);
        vert = [vert, temp];
    end
    
    % calculate the hausdorff-distance between the vertices and the point
    % cloud
    val = -inf;
    
    for i = 1:size(vert,2)
       temp = points - vert(:,i) * ones(1,size(points,2));
       len = sum(temp.^2,1);
       valTemp = min(len);
       val = max(val,sqrt(valTemp));
    end

% ------------------------------ END OF CODE ------------------------------
