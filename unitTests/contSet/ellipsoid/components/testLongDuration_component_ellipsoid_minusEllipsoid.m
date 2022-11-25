function res = testLongDuration_ellipsoid_orEllipsoidIA
% testLongDuration_ellipsoid_orEllipsoidIA - unit test function of testLongDuration_ellipsoid_orEllipsoidIA
%
% Syntax:  
%    res = testLongDuration_ellipsoid_orEllipsoidIA
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

% Author:       Victor Gassmann
% Written:      17-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        try
            % for outer approx, both need to be degenerate, otherwise "and"
            % will throw error 
            E1 = ellipsoid.generateRandom(i,false);
            E2 = ellipsoid.generateRandom(i,false);
            % test "direct syntax
            E1o = E1-E2;
            E1i = minus(E1,E2,'i');
            [U,S,~] = svd(E2.Q);
            m = ceil(dim(E2)*rand);
            S(m,m) = 0;
            % generate degenerate ellipsoid contained in E2
            E3 = ellipsoid(U*S*U',E2.q);
            E2i = minus(E1,E3,'i');
            if isempty(E1o) 
                if ~(E1==E2) && in(E1,E2)
                    res = false;
                    break;
                end
                continue;
            end
            if ~in(E1o,E1i) || ~in(E1o,E2i)
                res = false;
                break;
            end
            E3i = minus(E2,E3,'i');
            if isempty(E3i)
                res = false;
                break;
            end
            % test cell array syntax
            E2o = minus(E1,{E2});
            if ~(E1o==E2o)
                res = false;
                break;
            end
            % check points
            Y1 = randPoint(E1,nPoints,'extreme');
            Y2 = randPoint(E2,nPoints,'extreme');
            Y  = sumPoints(Y1,-Y2);
            %check if Y \in Eo
            if ~in(E1o,Y)
                res = false;
                break;
            end
        catch
            if strcmp(ME.identifier,'CORA:solverIssue')
                disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                continue;
            end
            rethrow(ME);
        end
    end
    if ~res
        break;
    end
end
%------------- END OF CODE --------------