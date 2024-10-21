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

% Authors:       Niklas Kochdumper
% Written:       03-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% this test runs forever
res = true; return
splits = 4;

% Random Tests --------------------------------------------------------

% define set representations that are tested
sets = {@conPolyZono.generateRandom,...
        @polyZonotope.generateRandom, ...
        @zonotope.generateRandom, ...
        @conZonotope.generateRandom, ...
        @ellipsoid.generateRandom, ...
        @capsule.generateRandom};

% loop over all test cases
for i = 1:2
    
    % generate random constrained polynomial zonotope
    cPZ1 = conPolyZono.generateRandom('Dimension',2,'NrIndGenerators',0);
    
    % loop over all other set representations
    for j = 1:length(sets)
        
        found = false;
        
        % loop unit two sets which really intersect are obtained
        for h = 1:3
        
            % generate random object of the current set representation
            if j == 1 || j == 2  % conPolyZono, polyZonotope
                S = sets{i}('Dimension',2,'NrIndGenerators',0);
            else
                S = sets{i}('Dimension',2);
            end
            
            cPZ2 = conPolyZono(S);

            % compute intersection
            cPZ = cPZ1 & S;

            % get random points
            try
                points = randPoint(cPZ,10,'extreme');
            catch
                continue 
            end
            
            found = true;
            break;
        end
        
        % check if two intersecting sets could be generated
        if ~found
            continue; 
        end
        
        % compute polytope enclosure of the two conPolyZono objects
        pgon1 = polygon(cPZ1,splits);
        pgon2 = polygon(cPZ2,splits);
        
        pgon = pgon1 & pgon2;
        
        % check if all points are inside the intersection of the
        % polygon enclosures
        assert(contains(pgon,points))
    end
end

% ------------------------------ END OF CODE ------------------------------
