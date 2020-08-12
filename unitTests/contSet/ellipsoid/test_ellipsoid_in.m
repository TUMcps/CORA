function res = test_ellipsoid_in
% test_ellipsoid_in - unit test function of in
%
% Syntax:  
%    res = test_ellipsoid_in
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
% Written:      15-October-2019
% Last update:  07-August-2020
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
condTOL = 1e4;
for i=10:5:15
    for j=1:nRuns
        while 1
            % make sure ellipsoid is not squished
            E1 = ellipsoid.generateRandom(false,i);
            if cond(E1.Q)<=condTOL
                break;
            end
        end
        E2 = ellipsoid(E1.Q*0.7,E1.q);
        if ~in(E1,E2)
            res = false;
            break;
        end
        while 1
            E1 = ellipsoid.generateRandom(false,i);
            E2 = ellipsoid.generateRandom(false,i);
            if cond(E1.Q)<=condTOL && cond(E2.Q)<=condTOL
                break;
            end
        end
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
    disp('test_ellipsoid_in successful');
else
    disp('test_ellipsoid_in failed');
end
%------------- END OF CODE --------------