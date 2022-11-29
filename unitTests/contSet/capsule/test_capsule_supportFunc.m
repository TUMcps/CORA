function res = test_capsule_supportFunc
% test_capsule_supportFunc - unit test function of supportFunc
%
% Syntax:  
%    res = test_capsule_supportFunc
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

% Analytical test: generator aligned on x_1, unit length; unit radius
res = true;

% empty set
C_e = capsule();
res = supportFunc(C_e,[1;1],'upper') == -Inf && supportFunc(C_e,[1;1],'lower') == Inf;

% loop over different dimensions
for n = 2:4:30
    
    % instantiate capsule
    c = zeros(n,1);
    g = [1;zeros(n-1,1)];
    r = 1;
    C = capsule(c,g,r);
    
    % support function should return
    %    -2|2 for first basis vector (lower|upper)
    %    -1|1 for all other basis vectors (lower|upper)
    
    % loop over basis vectors
    id = eye(n);
    for e=1:n
        % compute support function in direction of all basis vectors
        basisvector = id(:,e);
        [val_upper,vec_upper] = supportFunc(C,basisvector,'upper');
        [val_lower,vec_lower] = supportFunc(C,basisvector,'lower');
        
        % compare to analytical solution
        if e==1 && (val_lower ~= -2 || val_upper ~= 2 || ...
            	~all(vec_lower == [-2;zeros(n-1,1)]) || ~all(vec_upper == [2;zeros(n-1,1)]) )
            res = false; return
        elseif e~=1 && (val_lower ~= -1 || val_upper ~= 1 || ...
            	~all(vec_lower == -basisvector) || ~all(vec_upper == basisvector) )
            res = false; return
        end
    end

end

%------------- END OF CODE --------------