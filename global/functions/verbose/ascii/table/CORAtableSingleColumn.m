classdef CORAtableSingleColumn < CORAtable
% CORAtableSingleColumn - CORAtable with single column
%
% Syntax:
%    table = CORAtableSingleColumn(heading)
%    table = CORAtableSingleColumn(heading,design,varargin)
%
% Inputs:
%    heading - char
%    design - char
%    varargin - name-value pairs
%        <'ColumnWidths',colWidths> - numeric, column width
%
% Outputs:
%    table - CORAtableSingleColumn object
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
    function table = CORAtableSingleColumn(heading,varargin)
        % 1. parse input
        design = setDefaultValues({'modern'},varargin(1:min(1,numel(varargin))));
        [heading,design] = aux_parseInput(heading,design);

        % 2. call constructor
        table = table@CORAtable(design,{heading},{'s'},varargin{2:end});
    end

    function printContentRow(table,content,varargin)
        % build content
        content = sprintf(content,varargin{:});
        
        % make sure it fits cell width to avoid right alignment
        spacing = repmat(' ',max(0,table.colWidths(1)-numel(content)));
        contentAligned = sprintf('%s%s',content,spacing);

        % print row
        printContentRow@CORAtable(table,{contentAligned});   
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
