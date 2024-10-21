function res = test_CORAtableParameters
% test_CORAtableParameters - unit test function for CORAtableParameters
%
% Syntax:
%    res = test_CORAtableParameters
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       20-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% designs
designs = {'modern', 'minimalistic'};

% test each design
for i=1:numel(designs)
    design = designs{i};

    % init table
    table = CORAtableParameters('Parameter List',design);

    % print table
    fprintf('Design: %s\n',design)
    table.printHeader()
    table.printContentRow('Param1','Value1')
    table.printContentRow('Param2')
    table.printContentRow('Param2-1','Value2-1','s',2)
    table.printContentRow('Param2-1','Value2-2','s',2)
    table.printMidBoundaryRow();
    table.printContentRow('Param3','Value3')
    table.printFooter();
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
