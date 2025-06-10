function linsys_full = full(linsys)
% full - Converts system matrices/vectors to full representation
%
% Syntax:
% linsys_full = full(linsys)
%
% Inputs:
%    linsys - linearSys object
%
% Outputs:
%    linsys_full - linearSys object with full system matrices/vectors
%
% Example:
%    A = rand(10,10);
%    A(A < 0.8) = 0;
%    B = rand(10,10);
%    B(B < 0.8) = 0;
%    sys = linearSys(sparse(A),sparse(B));
%    sys_full = full(sys);
%    disp(sys);
%    disp(sys_full);
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

A = full(linsys.A);
B = full(linsys.B);
c = full(linsys.c);
C = full(linsys.C);
D = full(linsys.D);
k = full(linsys.k);
E = full(linsys.E);
F = full(linsys.F);

linsys_full = linearSys(linsys.name,A,B,c,C,D,k,E,F);

end

% ------------------------------ END OF CODE ------------------------------
