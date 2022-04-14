function res = testLongDuration_conZonotope_and
% testLongDuration_conZonotope_and - unit test function for intersection
%    of a constrained zonotope with other sets
%
% Syntax:  
%    res = testLongDuration_conZonotope_and
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

% TEST 1: conZonotope (random) --------------------------------------------

% loop over different dimensions
for j = 2:3
    
    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    poly1 = mptPolytope(V');

    % Generate random polytope vertices 1
    points = rand(j,100);
    ind = convhulln(points');
    ind = unique(ind(:,1),'stable');
    V = points(:,ind);
    poly2 = mptPolytope(V');

    % calculate constrained zonontope intersection
    cZono1 = conZonotope(poly1);
    cZono2 = conZonotope(poly2);
    zonoInt = cZono1 & cZono2;
    V = vertices(zonoInt);

    % calculate vertices from polytope interesection
    polyInt = poly1 & poly2;
    V_ = vertices(polyInt);

    % plot the result
%     if j == 2
%         plot(cZono1,[1,2],'r');
%         hold on
%         plot(cZono2,[1,2],'b');
%         plot(zonoInt,[1,2],'g','Filled',true,'EdgeColor','none');
%         plot(V(1,:),V(2,:),'.k','MarkerSize',12);
%     end

    % check correctness
    if ~isempty(V_)
        for i = 1:size(V_,2)
           if ~ismembertol(V_(:,i)',V',1e-10,'ByRows',true)
              file_name = strcat('testLongDuration_conZonotope_intersection_2_', ...
                                 datestr(now,'mm-dd-yyyy_HH-MM'));
                  
              file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                                   file_name);
                           
              save(file_path, 'poly1', 'poly2')
              error('Test 1 (conZonotope random) failed!'); 
           end
        end
    end
end

res = true;


%------------- END OF CODE --------------