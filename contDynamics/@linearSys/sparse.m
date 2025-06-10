function linsys_sparse = sparse(linsys)
% sparse - Converts system matrices/vectors to sparse representation
%
% Syntax:
% linsys_sparse = sparse(linsys)
%
% Inputs:
%    linsys - linearSys object
%
% Outputs:
%    linsys_sparse - linearSys object with sparse system matrices/vectors
%
% Example:
%    A = rand(10,10);
%    A(A < 0.8) = 0;
%    B = rand(10,10);
%    B(B < 0.8) = 0;
%    sys = linearSys(A,B);
%    sys_sparse = sparse(sys);
%    disp(sys);
%    disp(sys_sparse);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       02-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

A = sparse(linsys.A);
B = sparse(linsys.B);
c = sparse(linsys.c);
C = sparse(linsys.C);
D = sparse(linsys.D);
k = sparse(linsys.k);
E = sparse(linsys.E);
F = sparse(linsys.F);

linsys_sparse = linearSys(linsys.name,A,B,c,C,D,k,E,F);

end

% ------------------------------ END OF CODE ------------------------------
