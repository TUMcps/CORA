function res = testLongDuration_conZonotope_mptPolytope
% testLongDuration_conZonotope_mptPolytope - unit test function for
%    conversion between constrained zonotopes and polytopes
%
% Syntax:  
%    res = testLongDuration_conZonotope_mptPolytope
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% TEST 1: 2D --------------------------------------------------------------

for j = 1:5
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
    V1 = vertices(cZ);

    % Convert back to a mptPolytope
    P = mptPolytope(cZ);

    % Calculate vertices
    V2 = vertices(P);

    % plot the result
%     plot(cZono,[1,2],'FaceColor','b');
%     hold on
%     plot(poly,[1,2],'r');
%     plot(V(1,:),V(2,:),'.k','MarkerSize',12);


    % Check for correctness
    if ~compareMatrices(V,V1,1e-10)
       file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                          datestr(now,'mm-dd-yyyy_HH-MM'));
              
       file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                            file_name);
                       
       save(file_path, 'cZ')
       throw(CORAerror('CORA:testFailed'));
    end

    if ~compareMatrices(V,V2,1e-10)
       file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                          datestr(now,'mm-dd-yyyy_HH-MM'));
              
       file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                            file_name);
                       
       save(file_path, 'P')
       throw(CORAerror('CORA:testFailed'));
    end
end


% TEST 2: 3D --------------------------------------------------------------

for j = 1:5
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
    V1 = vertices(cZ);

    % Convert back to a mptPolytope
    P = mptPolytope(cZ);

    % Calculate vertices
    V2 = vertices(P);

    % Check for correctness
    if ~compareMatrices(V,V1,1e-10)
       file_name = strcat('testLongDuration_conZonotope_polytope_2_', ...
                          datestr(now,'mm-dd-yyyy_HH-MM'));
              
       file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                            file_name);
                       
       save(file_path, 'cZ')
       throw(CORAerror('CORA:testFailed'));
    end

    if ~compareMatrices(V,V2,1e-10)
       file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                          datestr(now,'mm-dd-yyyy_HH-MM'));
              
       file_path = fullfile(CORAROOT, 'unitTests', 'failedTests', ...
                            file_name);
                       
       save(file_path, 'P')
       throw(CORAerror('CORA:testFailed'));
    end
end

%------------- END OF CODE --------------