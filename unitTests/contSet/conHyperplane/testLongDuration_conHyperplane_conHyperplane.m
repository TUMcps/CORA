function res = testLongDuration_conHyperplane_conHyperplane
% testLongDuration_conHyperplane_conHyperplane - unit test function of
%    conHyperplane (constructor)
%
% Syntax:  
%    res = testLongDuration_conHyperplane_conHyperplane
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
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-12;

res = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(25);
    % number of constraints
    nrCon = randi(25);
    
    % random normal vector, offset, constraint matrix, constraint vector
    a = randn(n,1);
    b = randn(1);
    hs = halfspace(a,b);
    C = randn(nrCon,n);
    d = randn(nrCon,1);
    
    % admissible initializations
    % only halfspace
    hyp = conHyperplane(hs);
    if ~isequal(hyp.h,hs)
        res = false; break;
    end

    % normal vector and offset
    hyp = conHyperplane(a,b);
    if ~isequal(hyp.h,halfspace(a,b))
        res = false; break;
    end

    % halfspace and constraint matrix, constraint vector
    try
        hyp = conHyperplane(hs,C,d);
        if ~isequal(hyp.h,hs) || any(any(abs(hyp.C - C) > tol)) || ...
                any(abs(hyp.d - d) > tol)
            res = false; break;
        end
    catch ME
        if strcmp(ME.identifier,'CORA:wrongInputInConstructor')
            % find a new random example
            continue;
        else
            rethrow(ME);
        end
    end
        
    % normal vector, offset, and constraint matrix, constraint vector
    hyp = conHyperplane(a,b,C,d);
    if ~isequal(hyp.h,halfspace(a,b)) || ...
            any(any(abs(hyp.C - C) > tol)) || any(abs(hyp.d - d) > tol)
        res = false; break;
    end
    
    
    % wrong initializations
    a_plus1 = randn(n+1,1);
    b_vec = randn(n+1,1);
    hs_plus1 = halfspace(a_plus1,b);
    C_plus1 = randn(nrCon,n+1);
    d_plus1 = randn(nrCon+1,1);
    
    % offset as vector
    try
        hyp = conHyperplane(a,b_vec); % <- should throw error here
        res = false; break;
    end
    
    % C and d do not fit halfspace
    try
        hyp = conHyperplane(hs_plus1,C,d); % <- should throw error here
        res = false; break;
    end
    try
        hyp = conHyperplane(hs,C_plus1,d); % <- should throw error here
        res = false; break;
    end 
    
    % C does not fit d
    try
        hyp = conHyperplane(hs,C,d_plus1); % <- should throw error here
        res = false; break;
    end 
    
    % a does not fit C
    try
        hyp = conHyperplane(a,b,C_plus1,d); % <- should throw error here
        res = false; break;
    end 
    
    % a does not fit d
    try
        hyp = conHyperplane(a,b,C,d_plus1); % <- should throw error here
        res = false; break;
    end 
    
    % too many input arguments
    try
        hyp = conHyperplane(a,b,C,d,d); % <- should throw error here
        res = false; break;
    end 
end


if res
    disp('testLongDuration_conHyperplane successful');
else
    disp('testLongDuration_conHyperplane failed');
end

%------------- END OF CODE --------------