function ls = levelSet(E)
% levelSet - Converts an ellipsoid to a level set
%
% Syntax:
%    ls = levelSet(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    ls - levelSet object
%
% Example: 
%    E = ellipsoid([4 2;2 4],[1;1]);
%    ls = levelSet(E);
% 
%    figure; hold on; xlim([-2,4]); ylim([-2,4]);
%    plot(ls,[1,2],'r');
%    plot(E,[1,2],'EdgeColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Niklas Kochdumper
% Written:       09-April-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'}});

% construct level set
x = sym('x',[size(E.Q,1),1]);
eq = transpose(x - E.q) * inv(E.Q) * (x - E.q) - 1;

ls = levelSet(eq,x,'<=');
    
% ------------------------------ END OF CODE ------------------------------
