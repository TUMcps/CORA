function res = test_zonotope_enclose
% test_enclose - unit test function of enclose
%
% Syntax:  
%    res = test_zonotope_enclose
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical Test ------------------------------------------------------

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9, 10, 11; 12, 13, 14]);

% obtain enclosing zonotope
Z3 = enclose(Z1,Z2);

% obtain zonotope matrix
Zmat = Z3.Z;

% true result
true_mat = [5, 6, 7, -4, -4, -4, 4; ...
            8.5, 9.5, 10.5, -3.5, -3.5, -3.5, 8];

% check result
res_val = all(all(Zmat == true_mat));


% 2. Random Tests ---------------------------------------------------------

dims = 2:4;
testsPerDim = 50;
ptsPerLine = 10;

% box has to be the same as conversion to interval
for d=1:length(dims)
    for test=1:testsPerDim
        % create two random zonotopes
        nrOfGens = randi([5,15],1,1);
        c1 = -20*ones(dims(d),1);
        G1 = -1+2*rand(dims(d),nrOfGens);
        Z1 = zonotope(c1,G1);
        c2 = 20*ones(dims(d),1);
        G2 = -1+2*rand(dims(d),nrOfGens);
        Z2 = zonotope(c2,G2);
        
        % compute enclosure
        Zenc = enclose(Z1,Z2);

        % random points in Z1 or Z2
        p1 = randPoint(Z1);
        p2 = randPoint(Z2);
        % connect points by line, all have to be in enclosure
        pts = p1 + (p2-p1) .* linspace(0,1,ptsPerLine);

        % random points have to be also in Zenc
        res_rand(d,test) = all(containsPoint(Zenc,pts));
    end
end


% add results
res = res_val && all(all(res_rand));

if res
    disp('test_enclose successful');
else
    disp('test_enclose failed');
end

%------------- END OF CODE --------------
