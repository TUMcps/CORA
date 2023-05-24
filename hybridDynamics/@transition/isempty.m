function res = isempty(trans)
% isempty - checks if a transition object is empty
%
% Syntax:  
%    res = isequal(trans)
%
% Inputs:
%    trans - transition object
%
% Outputs:
%    res - true/false
%
% Example: 
%    % guard set
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
%
%    % transition
%    trans = transition(guard,reset,1);
%
%    % comparison
%    res = isempty(trans)
%    res = isempty(transition())
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init logical array
[r,c] = size(trans);
res = false(r,c);

% loop over all objects
for i=1:r
    for j=1:c
        % check target (has to be given)
        targ = trans(i,j).target;
        if isnumeric(targ) && isempty(targ)
            res(i,j) = true;
        end
    end
end

%------------- END OF CODE --------------