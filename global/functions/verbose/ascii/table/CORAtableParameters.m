classdef CORAtableParameters < CORAtableSingleColumn
% CORAtableParameters - CORAtable listing parameters
%
% Syntax:
%    table = CORAtableParameters(heading)
%    table = CORAtableParameters(heading,design,varargin)
%
% Inputs:
%    heading - char
%    design - char
%    varargin - name-value pairs
%        <'ColumnWidths',colWidths> - numeric, column width
%
% Outputs:
%    table - CORAtableParameters object
%
% Example: 
%    % init table
%    table = CORAtableParameters('Parameter List');
%
%    % print table
%    table.printHeader()
%    table.printContentRow('Param1','Value1')
%    table.printContentRow('Param2')
%    table.printContentRow('Param2-1','Value2-1','s',2)
%    table.printContentRow('Param2-1','Value2-2','s',2)
%    table.printContentRow('Param3','Value3')
%    table.printFooter();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAtable

% Authors:       Tobias Ladner
% Written:       19-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods
    function table = CORAtableParameters(heading,varargin)
        table = table@CORAtableSingleColumn(heading,varargin{:});
    end

    function printContentRow(table,name,varargin)
        % parse input
        narginchk(2,5)
        [value,format,level] = setDefaultValues({'','s',1},varargin);

        % obtain format
        
        % computed intend
        indent = repmat('  ',1,level-1);
        
        printContentRow@CORAtableSingleColumn(table,['%s- %s: %' format],indent,name,value);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [heading,design] = aux_parseInput(heading,design)

    % check input args
    inputArgsCheck({ ...
        {'heading','att','char'};
        {'design','att','char'}; % is checked in CORAtable
    })
end

% ------------------------------ END OF CODE ------------------------------
