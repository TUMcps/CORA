function res = testLongDuration_capsule_supportFunc
% testLongDuration_capsule_supportFunc - unit test function of supportFunc
%
% Syntax:  
%    res = testLongDuration_capsule_supportFunc
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

% Author:       Mark Wetzlinger
% Written:      11-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
            break
        end
    end
    
    if ~res
        break;
    end
    
end


if res
    disp('test_supportFunc successful');
else
    disp('test_supportFunc failed');
end

%------------- END OF CODE --------------