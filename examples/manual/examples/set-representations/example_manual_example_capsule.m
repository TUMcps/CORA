function example_manual_example_capsule()
% example_manual_example_capsule - example from the manual demonstrating 
% the capsule example from the manual
%
% Syntax:
%   example_manual_example_capsule()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct a capsule
c = [1;2];
g = [2;1];
r = 1;

C1 = capsule(c,g,r)

% linear map of a capsule
A = [0.5 0.2; -0.1 0.4];
C2 = A*C1;

% shift the center of a capsule
s = [0;1];
C3 = C2 + s;

% check capsule-in-capsule containment
res1 = contains(C1,C2);
res2 = contains(C1,C3);

disp(['C2 in C1?: ',num2str(res1)]);
disp(['C3 in C1?: ',num2str(res2)]);

% visualization
figure; hold on
plot(C1,[1,2],'r');
plot(C2,[1,2],'g');
plot(C3,[1,2],'b');

% plot --------------------------------------------------------------------

enlargeAxis(1.2)
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
