function res = isequal(trans1,trans2)
% isequal - checks if two transitions are equal by comparing the guard
%    sets, reset functions, target locations, and synchronization labels
%
% Syntax:  
%    res = isequal(trans1,trans2)
%
% Inputs:
%    trans1 - transition object
%    trans2 - transition object
%
% Outputs:
%    res - true/false
%
% Example: 
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);
%
%    % transition
%    trans1 = transition(guard,reset1,1);
%    trans2 = transition(guard,reset2,1);
%
%    % comparison
%    res = isequal(trans1,trans2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% easy checks first to facilitate quick exits
res = true;

% target location
if any(size(trans1.target) ~= size(trans2.target)) ...
        || ~all(trans1.target == trans2.target)
    res = false; return
end

% synchronization label
if ~strcmp(trans1.syncLabel,trans2.syncLabel)
    res = false; return
end

% reset function (conversion from nonlinear to linear not checked)

% check fields (either both linear or both nonlinear)
reset1fields = fields(trans1.reset);
reset2fields = fields(trans2.reset);
if length(reset1fields) ~= length(reset2fields)
    res = false; return
end
idx2 = zeros(length(reset1fields),1);
for i=1:length(reset1fields)
    temp = cellfun(@(x) strcmp(x,reset1fields{i}),reset2fields,'UniformOutput',true);
    if isempty(temp)
        % field from reset of trans1 does not exist in reset of trans2
        res = false; return
    end
    % save idx for value comparison
    idx2(i) = find(temp == true,1,'first');
end
% established that both reset functions have the same fields
for i=1:length(reset1fields)
    fieldname = reset1fields{i};
    switch fieldname
        case 'A'
            % state mapping
            if any(size(trans1.reset.A) ~= size(trans2.reset.A)) ...
                    || ~all(all(trans1.reset.A == trans2.reset.A))
                res = false; return
            end
        case 'B'
            % input mapping
            if any(size(trans1.reset.B) ~= size(trans2.reset.B)) ...
                    || ~all(all(trans1.reset.B == trans2.reset.B))
                res = false; return
            end
        case 'c'
            % constant offset
            if any(size(trans1.reset.c) ~= size(trans2.reset.c)) ...
                    || ~all(all(trans1.reset.c == trans2.reset.c))
                res = false; return
            end
        case 'f'
            % function handle for nonlinear reset function
            if ~isequal(trans1.reset.f,trans2.reset.f)
                res = false; return
            end
        case 'stateDim'
            % dimension of state
            if trans1.reset.stateDim ~= trans2.reset.stateDim
                res = false; return
            end
        case 'inputDim'
            if trans1.reset.inputDim ~= trans2.reset.inputDim
                res = false; return
            end
        case 'hasInput'
            if trans1.reset.hasInput ~= trans2.reset.hasInput
                res = false; return
            end
    end
end

% guard set
if ~isequal(trans1.guard,trans2.guard)
    res = false; return
end

%------------- END OF CODE --------------