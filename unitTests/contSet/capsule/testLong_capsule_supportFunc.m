function res = testLong_capsule_supportFunc
% testLong_capsule_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = testLong_capsule_supportFunc
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

% Authors:       Mark Wetzlinger
% Written:       11-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random test: random generator in unit sphere, unit radius
res = true;

for n=2:4:30
    
    % instantiate capsule
    c = zeros(n,1);
    g = randomPointOnSphere(n);
    r = 1;
    C = capsule(c,g,r);
    
    % directions for support function
    numDirs = 1000;
    sections = nthroot(numDirs,n-1);
    dirs = randEqdistDirections(n,sections);
    
    for d=1:size(dirs,2)
        val = supportFunc(C,dirs(:,d));
        if abs(val) > 2
            res = false;
            return
        end
    end    
end

% ------------------------------ END OF CODE ------------------------------
