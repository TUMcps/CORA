function res = testLong_conPolyZono_reduceConstraints
% testLong_conPolyZono_reduceConstraints - unit test function for
%    the constraint reduction of constrained polynomial zonotopes
%
% Syntax:  
%    res = testLong_conPolyZono_reduceConstraints()
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
% See also: conPolyZono/reduceConstraints

% Author:       Niklas Kochdumper
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
splits = 4;

% loop over all test cases
for i = 1:10
    
    % generate random constrained polynomial zonotopes with constraints
    con = 0;
    while con < 1
        cPZ1 = conPolyZono.generateRandom('Dimension',2);
        con = size(cPZ1.A,1);
    end

    % reduce constraints
    nrCon = randi([0,con-1]);
    cPZ = reduceConstraints(cPZ1,nrCon);

    % compute random points inside the original set
    points = randPoint(cPZ1,10,'extreme');
    
    % check if all points are inside polygon enclosures
    pgon = polygon(cPZ,splits);

    if ~contains(pgon,points)
        throw(CORAerror('CORA:testFailed'));
    end
end

%------------- END OF CODE --------------