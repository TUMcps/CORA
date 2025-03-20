function res = testLong_conPolyZono_and
% testLong_conPolyZono_and - unit test function for the
%    intersection of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_and()
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
% See also: conPolyZono/and

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       03-February-2021
% Last update:   05-March-2025 (TL, rewrote)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random Tests --------------------------------------------------------

% define set representations that are tested
sets = {@conPolyZono.generateRandom,...
        @polyZonotope.generateRandom, ...
        @zonotope.generateRandom, ...
        @conZonotope.generateRandom, ...
        @ellipsoid.generateRandom, ...
        @capsule.generateRandom};
numTestsPerSet = 5;

% loop over all other set representations
for i = 1:numel(sets)
    % generate multiple test cases per set
    for j = 1:numTestsPerSet
        fprintf('i=%i/%i (%s), j=%i/%i..\n',i,numel(sets),func2str(sets{i}),j,numTestsPerSet)
    
        % generate random conPolyZono for which a point can be determined
        p1 = [];
        while isempty(p1)
            % generate random constrained polynomial zonotope
            cPZ1 = conPolyZono.generateRandom('Dimension',2,'NrIndGenerators',0);
            
            % sample point
            p1 = randPoint(cPZ1);
        end
    
        % generate second random set for which a point can be determined
        p2 = [];
        while isempty(p2)
            % generate random constrained polynomial zonotope
            if i <= 2  % conPolyZono, polyZonotope
                S2 = sets{i}('Dimension',2,'NrIndGenerators',0);
            else
                S2 = sets{i}('Dimension',2);
            end
            
            % sample point
            p2 = randPoint(S2);
        end
    
        % shift S2 such that points overlap -> p1 must be in intersection
    
        % compute distance between determined points
        dist = p1-p2;
    
        % shift S2
        S2 = S2 + dist;
    
        % compute intersection (p1 is in both sets)
        cPZ_inter = cPZ1 & S2;
        
        % compute polygon enclosure
        splits = 4;
        pgon_inter = polygon(cPZ_inter,splits);
        
        % check if p1 is contained in the intersection
        assertLoop(contains(pgon_inter,p1),i,j)
    end
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
