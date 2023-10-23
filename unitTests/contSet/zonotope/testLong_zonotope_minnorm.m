function res = testLong_zonotope_minnorm
% testLong_zonotope_minnorm - unit test function of minnorm
%
% Syntax:
%    res = testLong_zonotope_minnorm
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

% Authors:       Victor Gassmann
% Written:       15-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% loop over dimension
for i=2:7

    % loop over number of generators
    for j=i:5:20

        % init random zonotope
        Z = zonotope([zeros(i,1),randn(i,j)]);

        % compute minnorm
        val = minnorm(Z);

        % evaluate support function for random unit directions
        L = eq_point_set(i-1,2*j);
        for k=1:length(L)
            sF = supportFunc(Z,L(:,k));
            % ensure that val <= suppfnc(Z,l)
            if val > sF && ~withinTol(val,sF)
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
