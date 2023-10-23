function completed = example_capsule()
% example_capsule - example demonstrating set based computation with
%                   capsules
%
% Syntax:
%    completed = example_capsule()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ---
% Written:       ---
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
C2 = A * C1;

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
plot(C1,[1,2]);
plot(C2,[1,2],'Color',colorblind('r'));
plot(C3,[1,2],'k');

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
