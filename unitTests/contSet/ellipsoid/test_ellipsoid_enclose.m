function res = test_ellipsoid_enclose
% test_ellipsoid_enclose - unit test function of enclose
%
% Syntax:  
%    res = test_ellipsoid_enclose
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

% Author:       Victor Gaﬂmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%We need at least one non-degenerate ellipsoid to ensure
%the union of E1,E2 is not the empty set.
nRuns = 2;
res = true;
for i=10:5:20
    for j=1:nRuns
        E1 = ellipsoid.generateRandom(false,i);
        %although enclose supports one degenerate ellipsoid, 'in' does not at the
        %moment, thus we only test full-dimensional ellipsoids
        %E2 = ellipsoid.generateRandom(n);
        E2 = ellipsoid.generateRandom(false,i);
        E = enclose(E1,E2);
        if ~in(E,E1) || ~in(E,E2)
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end
if res
    disp('test_ellipsoid_enclose successful');
else
    disp('test_ellipsoid_enclose failed');
end


%------------- END OF CODE --------------
