function res = test_CORAtableSummary
% test_CORAtableSummary - unit test function for CORAtable
%
% Syntax:
%    res = test_CORAtableSummary
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

% Authors:       Lukas Koller
% Written:       01-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% designs
designs = {'single', 'double', 'modern', 'minimalistic', 'ascii','latex'};

% test each design
for i=1:numel(designs)
    design = designs{i};
    fprintf('Design: %s\n\n',design)

    % Table 1 ---

    % init table
    table = CORAtable(design,{'Heading 1','Heading 2'}, ...
        {'s','sum{%.3e & %.3e}'});

    % print table
    table.printHeader()
    table.printContentRow({'v11',0.0:0.01:0.4});
    table.printContentRow({'v21',0.0:0.02:0.5});
    table.printMidBoundaryRow();
    table.printContentRow({'v31',0.0:0.03:0.6});
    table.printFooter();

end

% ---

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
