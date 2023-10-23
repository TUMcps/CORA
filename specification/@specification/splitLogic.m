function [specNonLogic,specLogic] = splitLogic(spec)
% splitLogic - split into temporal logic and non-logic specifications
%
% Syntax:
%    [spec,specLogic] = splitLogic(spec)
%
% Inputs:
%    spec - specification object
%
% Outputs:
%    specNonLogic - non-temporal logic specifications
%    specLogic - temporal logic specifications
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper
% Written:       17-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

specNonLogic = []; specLogic = [];

for i = 1:size(spec,1)
    if strcmp(spec(i,1).type,'logic')
        % temporal logic specifications
        specLogic = add(specLogic,spec(i,1));
    else
        % non-temporal logic specifications
        specNonLogic = add(specNonLogic,spec(i,1));
    end
end

% ------------------------------ END OF CODE ------------------------------
