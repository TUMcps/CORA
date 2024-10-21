function SpS = spectraShadow(I)
% spectraShadow - Converts an interval to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    SpS - spectraShadow object
%
% Example:
%    I = interval([-3.14;-4.20],[2.718;0.69]);
%    SpS = spectraShadow(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if not a matrix set
n = dim(I);
if numel(n) > 1
    throw(CORAerror('CORA:wrongValue','first','Interval must not be an n-d array with n > 1.'))
end

lower = I.infimum;
upper = I.supremum;


A0_diagonals = cell([1 n]);
Ai = cell([1 n]);
for i=1:n
    D_A0 = [upper(i) 0;0 -lower(i)];
    
    % We need several case-by-case scenarios for the cases where one of the
    % bounds is, well, unbounded
    if upper(i)==Inf
        D_A0(1,1) = 1;
    elseif upper(i)==-Inf
        D_A0(1,1) = -1;
    end
    if lower(i)==Inf
        D_A0(2,2) = -1;
    elseif lower(i)==-Inf
        D_A0(2,2) = 1;
    end
    A0_diagonals{i} = sparse(D_A0);
    
    D_Ai = [-1 0;0 1];
    % Do the same for the Ai matrices
    if upper(i)==Inf
        D_Ai(1,1) = 0;
    elseif upper(i)==-Inf
        D_Ai(1,1) = 0;
    end
    if lower(i)==Inf
        D_Ai(2,2) = 0;
    elseif lower(i)==-Inf
        D_Ai(2,2) = 0;
    end
    
    Ai{i} = blkdiag(sparse(2*(i-1),2*(i-1)),sparse(D_Ai),sparse(2*(n-i),2*(n-i)));
end
A0 = blkdiag(A0_diagonals{:});

SpS = spectraShadow([A0 cat(2,Ai{:})]);

% Additional properties
SpS.bounded.val = isBounded(I);
SpS.emptySet.val = representsa_(I,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(I);
SpS.center.val = center(I);

% ------------------------------ END OF CODE ------------------------------
