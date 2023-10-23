function res = testLong_polyZonotope_polyMap
% testLong_polyZonotope_polyMap - unit test function of polyMap
%
% Syntax:
%    res = testLong_polyZonotope_polyMap
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

% Authors:       Niklas Kochdumper
% Written:       23-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    for i = 1:100

        % generate random polynomial zonotope
        pZ = polyZonotope.generateRandom("NrGenerators",3);
        
        % generate random polynomial map
        map = polyZonotope.generateRandom("NrFactors",dim(pZ));
        
        coeff = map.G;
        E = map.E;
    
        % compute polynomial map
        pZres = polyMap(pZ,coeff,E);
    
        % test result for random points
        fac = size(pZ.E,1);
    
        for j = 1:100
                
            % determine random point the polynoimal zonotope
            alpha = randPoint(interval(-ones(fac,1), ones(fac,1)));
            p = pZ.c + sum(pZ.G.*prod(alpha.^pZ.E,1),2);
    
            % compare exact result with result obtained from polynomial
            % zonotope computation
            p1 = pZres.c + sum(pZres.G.*prod(alpha.^pZres.E,1),2);
            p2 = sum(coeff.*prod(p.^E,1),2);
    
            if norm(p1-p2) > 1e-5 && norm(p1-p2)/norm(p1) > 1e-5
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end

    res = true;
end

% ------------------------------ END OF CODE ------------------------------
