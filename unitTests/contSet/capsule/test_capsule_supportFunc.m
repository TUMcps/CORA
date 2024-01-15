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

% Analytical test: generator aligned on x_1, unit length; unit radius
res = true(0);

% empty capsule
C_e = capsule.empty(2);
res(end+1,1) = supportFunc(C_e,[1;1],'upper') == -Inf ...
    && supportFunc(C_e,[1;1],'lower') == Inf;

% 2D capsule without radius
C = capsule([1;-1],[0;0],3);
dir = [1;-2];
[valC,xC] = supportFunc(C,dir);
% check result
res(end+1,1) = withinTol(dir'*xC,valC) && contains(C,xC);


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
        [vals,vecs] = supportFunc(C,basisvector,'range');
        
        % compare to analytical solution
        if e==1 && ( ~withinTol(val_lower,-2) || ~withinTol(val_upper,2) ...
                || ~isequal(vals,interval(-2,2)) ...
                || ~compareMatrices([vec_lower,vec_upper],[[-2;zeros(n-1,1)],[2;zeros(n-1,1)]]) ...
                || ~compareMatrices(vecs,[[-2;zeros(n-1,1)],[2;zeros(n-1,1)]]) )
            throw(CORAerror('CORA:testFailed'));
        elseif e~=1 && (~withinTol(val_lower,-1) || ~withinTol(val_upper,1) ...
                || ~isequal(vals,interval(-1,1)) ...
                || ~compareMatrices([vec_lower,vec_upper],[-basisvector,basisvector]) ...
                || ~compareMatrices(vecs,[-basisvector,basisvector]) )
            throw(CORAerror('CORA:testFailed'));
        end
    end

end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
