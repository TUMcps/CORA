function res = testLong_conZonotope_reduce
% testLong_conZonotope_reduce - unit test function for order
%    reduction of constrained zonotopes
%
% Syntax:  
%    res = testLong_conZonotope_reduce
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      11-May-2018
% Last update:  05-December-2020
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% TEST 1: 2D --------------------------------------------------------------

methods = {'girard','combastel'};

for j = 1:length(methods)
    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a mptPolytope object from the vertices
    P = mptPolytope(V');

    % Convert to constrained zonotope object
    cZ = conZonotope(P);

    % Calculate vertices
    V = vertices(cZ);
    
    % Reduce all constraints for which the elimination does not result in
    % an over-approximation
    cZ = reduceConstraints(cZ);

    % Reduce the constrained zonotope
    nc = size(cZ.A,1);
    ng = size(cZ.Z,2)-1;
    n = size(cZ.Z,1);
    o = floor((ng-nc)/n)+1;
    
    cRed{1} = reduceConstraints(cZ,nc-1);   % reduce 1 constraint
    cRed{2} = reduceConstraints(cZ,nc-2);   % reduce 2 constraints
    cRed{3} = reduce(cZ,methods{j},o-1);    % reduce 1 generator
    cRed{4} = reduce(cZ,methods{j},o-4);    % reduce 2 generators
    
    % Check if the vertices of the original constrained zonotope are
    % located inside the reduced constrained zonotope
    for i = 1:length(cRed)
       
        % convert to polytope (for easy checks if point is inside set)
        P = mptPolytope(cRed{i});
        P = get(P,'P');
        
%         % plot the result
%         hold on
%         plot(cZ,[1,2],'r');
%         plot(cRed{i},[1,2],'b');
        
        % check if all vertices are located inside the set
        temp = P.A*V - P.b*ones(1,size(V,2));
        if ~all(all( temp < 0 | withinTol(temp,1e-10) ))
            throw(CORAerror('CORA:testFailed'));
        end
    end
end



% TEST 2: 3D --------------------------------------------------------------

methods = {'girard','combastel'};

for j = 1:length(methods)
    % Generate random polytope vertices
    points = rand(3,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a mptPolytope object from the vertices
    P = mptPolytope(V');

    % Convert to constrained zonotope object
    cZ = conZonotope(P);

    % Calculate vertices
    V = vertices(cZ);
    
    % Reduce all constraints for which the elimination does not result in
    % an over-approximation
    cZ = reduceConstraints(cZ);

    % Reduce the constrained zonotope
    nc = size(cZ.A,1);
    ng = size(cZ.Z,2)-1;
    n = size(cZ.Z,1);
    o = floor((ng-nc)/n)+1;
    
    cRed{1} = reduceConstraints(cZ,nc-1);    % reduce 1 constraint
    cRed{2} = reduceConstraints(cZ,nc-2);    % reduce 2 constraints
    cRed{3} = reduce(cZ,methods{j},o-1);     % reduce 1 generator
    cRed{4} = reduce(cZ,methods{j},o-2);     % reduce 2 generators
    
    % Chech if the vertices of the original constrained zonotope are
    % located inside the reduced constrained zonotope
    for i = 1:length(cRed)
       
        % convert to polytope (for easy checks if point is inside set)
        P = mptPolytope(cRed{i});
        P = get(P,'P');
        
        % check if all vertices are located inside the set
        temp = P.A*V - P.b*ones(1,size(V,2));
        if ~all(all( temp < 0 | withinTol(temp,1e-10) ))
            throw(CORAerror('CORA:testFailed'));
        end
    end
end


% Test 3: Redundant Constraints (analytical test) -------------------------

load('conZonotope_reduce.mat');

% remove all redundant constraints
cZred = reduceConstraints(cZ);

% calculate vertices
V = vertices(cZ)';
V_ = vertices(cZred)';

% compare with the real vertices
if ~compareMatrices(V,V_)
    res = false;
end

%------------- END OF CODE --------------