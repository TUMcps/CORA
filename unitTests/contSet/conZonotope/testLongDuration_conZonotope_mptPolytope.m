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
%    res - boolean 
%
% Example: 
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

res = false;

% TEST 1: 2D --------------------------------------------------------------

for j = 1:5
    % Generate random polytope vertices
    points = rand(2,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);

    % Construct a mptPolytope object from the vertices
    poly = mptPolytope(V');

    % Convert to constrained zonotope object
    cZono = conZonotope(poly);

    % Calculate vertices
    V1 = vertices(cZono);

    % Convert back to a mptPolytope
    poly = mptPolytope(cZono);

    % Calculate vertices
    V2 = vertices(poly);

    % plot the result
%     plot(cZono,[1,2],'b','Filled',true,'EdgeColor','none');
%     hold on
%     plot(poly,[1,2],'r');
%     plot(V(1,:),V(2,:),'.k','MarkerSize',12);


    % Check for correctness
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V1',1e-10,'ByRows',true)
           file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'cZono')
           error('Test 1 failed! Wrong conZonotope vertices!'); 
       end
    end

    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V2',1e-10,'ByRows',true)
           file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'poly')
          error('Test 1 failed! Wrong mptPolytope vertices!'); 
       end
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
    poly = mptPolytope(V');

    % Convert to constrained zonotope object
    cZono = conZonotope(poly);

    % Calculate vertices
    V1 = vertices(cZono);

    % Convert back to a mptPolytope
    poly = mptPolytope(cZono);

    % Calculate vertices
    V2 = vertices(poly);

    % Check for correctness
    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V1',1e-10,'ByRows',true)
           file_name = strcat('testLongDuration_conZonotope_polytope_2_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'cZono')
           error('Test 2 failed! Wrong conZonotope vertices!'); 
       end
    end

    for i = 1:size(V,2)
       if ~ismembertol(V(:,i)',V2',1e-10,'ByRows',true)
           file_name = strcat('testLongDuration_conZonotope_polytope_1_', ...
                              datestr(now,'mm-dd-yyyy_HH-MM'));
                  
           file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                                file_name);
                           
           save(file_path, 'poly')
           error('Test 2 failed! Wrong mptPolytope vertices!'); 
       end
    end
end




res = true;


%------------- END OF CODE --------------