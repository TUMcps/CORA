function res = test_CORAtable
% test_CORAtable - unit test function for CORAtable
%
% Syntax:
%    res = test_CORAtable
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
% Written:       19-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% designs
designs = {'single', 'double', 'modern', 'minimalistic', 'ascii', 'latex', 'html','markdown'};

% test each design
for i=1:numel(designs)
    design = designs{i};
    fprintf('Design: %s\n\n',design)

    % Table 1 ---

    % init table
    table = CORAtable(design,{'Heading 1','Heading 2'},{'s','s'});

    % print table
    table.printHeader()
    table.printContentRow({'v11','v12'});
    table.printContentRow({'v21','v22'});
    table.printMidBoundaryRow();
    table.printContentRow({'v31','v32'});
    table.printFooter();

    % Table 2 ---

    % init table
    table = CORAtable(design,{'i','Value'},{'i','.4e'},'ColumnWidths',[5,10]);

    % print table
    table.printHeader()
    table.printContentRow({1,1e-2});
    table.printContentRow({2,5e-3});
    table.printContentRow({3,2e-3});
    table.printFooter();
end

% ---

% formats should be without leading %
assertThrowsAs(@CORAtable,'CORA:wrongValue',"single",{'Heading 1'}, {'%s'});

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
