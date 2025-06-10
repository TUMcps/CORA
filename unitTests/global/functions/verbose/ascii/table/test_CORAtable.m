function res = test_CORAtable(varargin)
% test_CORAtable - unit test function for CORAtable
%
% Syntax:
%    res = test_CORAtable
%    res = test_CORAtable(designs)
%
% Inputs:
%    designs - cell array specifying designs to print
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

% parse input (only view specific designs)
alldesigns = {'single', 'double', 'modern', 'minimalistic', 'ascii', 'latex', 'html','markdown','csv'};
if nargin > 0
    designs = varargin{1};
    if strcmp(designs{1},'all')
        designs = alldesigns;
    end
else
    designs = alldesigns;
end

fprintf('Printing some tables in all designs: ---\n\n')

% test each design
for i=1:numel(designs)
    design = designs{i};
    fprintf('Design: %s --------------------------------------\n\n',design)

    % Table 1 ---
    disp('Table 1:')

    % init table
    coratable = CORAtable(design,{'#','Heading 1','Heading 2'},{'rownr','s','s'});

    % print table
    coratable.printHeader()
    coratable.printContentRow({[],'v11','v12'});
    coratable.printContentRow({[],'v21','v22'});
    coratable.printMidBoundaryRow();
    coratable.printContentRow({[],'v31','v32'});
    coratable.printFooter();

    % print to file
    filename = 'mytable1.txt';
    fprintf('Table 1: (printed to file: <a href="matlab: open %s">%s</a>)\n\n',filename,filename)

    % init table
    coratable = CORAtable(design,{'#','Heading 1','Heading 2'},{'rownr','s','s'},'FileName',filename);

    % print table
    coratable.printHeader()
    coratable.printContentRow({[],'v11','v12'});
    coratable.printContentRow({[],'v21','v22'});
    coratable.printMidBoundaryRow();
    coratable.printContentRow({[],'v31','v32'});
    coratable.printFooter();
    disp(' ')

    % Table 2 ---
    disp('Table 2:')

    % init table
    coratable = CORAtable(design,{'Time','i','Value','Detailed Statistics'},{'time','i','.4e','sum{%.4e & %.4e}'},'ColumnWidths',[8,5,10,20],'SaveContent',true);

    % print table
    coratable.printHeader()
    coratable.printContentRow({[],10,1e-2,[1.5e-2,1.3e-2,1.1e-2,1.2e-2]});
    coratable.printContentRow({[],7,5e-3,[2.1e-2,2.3e-2,2.7e-2,2.5e-2]});
    coratable.printContentRow({[],3,2e-3,[3.9e-2,4.0e-2,4.1e-2,3.7e-2]});
    coratable.printFooter();

    % reprint
    disp('Table 2 (reprint):')
    coratable.reprint(design)

    % reprint to file
    filename = 'mytable2.txt';
    fprintf('Table 2: (printed to file: <a href="matlab: open %s">%s</a>)\n\n',filename,filename)
    coratable.reprint(design,'FileName',filename)

    % Table 3 (from Matlab table) ---
    disp('Table 3:')

    LastName = ["Sanchez";"Johnson";"Zhang";"Diaz";"Brown"];
    Age = [38;43;38;40;49];
    Smoker = [true;false;true;false;true];
    Height = [71;69;64;67;64];
    Weight = [176;163;131;133;119];
    patients = table(LastName,Age,Smoker,Height,Weight);
    CORAtable.printTable(patients,design)
end

disp('-------------------------------------------------------')
fprintf('\n^ Scroll up to see tables (%i designs). ^\n\n',numel(designs))
fprintf('View designs: %s.\n', strjoin(cellfun(@(design) sprintf('<a href="matlab:test_CORAtable({''%s''})">%s</a>',design,design), [{'all'}, alldesigns], 'UniformOutput', false),', '))

% ---

% formats should be without leading %
assertThrowsAs(@CORAtable,'CORA:wrongValue',"single",{'Heading 1'}, {'%s'});

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
