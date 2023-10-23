function dispInput(varname)
% dispInput - Displays 'varname =' at the start of display() if not called
%    from other display() operation
%
% Syntax:
%    dispInput(varname)
%
% Inputs:
%    varname - name of variable (usually inputname(1))
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if called from display of hybridAutomaton
st = dbstack("-completenames");
callFromOtherDisplay = false;
if length(st) >= 3 && contains(st(3).file,'display')
    callFromOtherDisplay = true;
end

if ~callFromOtherDisplay
    fprintf(newline);
    disp([varname, ' =']);
    fprintf(newline);
end

% ------------------------------ END OF CODE ------------------------------
