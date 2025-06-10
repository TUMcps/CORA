classdef CORAtable < ASCIItable
% CORAtable - creates a table in the well-beloved  CORA format
%
% Syntax:
%    table = CORAtable(design,headings,formats,varargin)
%
% Inputs:
%    design - char, 'single', 'double', 'modern', 'minimalistic', 
%       'ascii', 'latex', 'html', 'markdown'
%    hvalues - cell array containing headings
%    formats - cell array containing cell formats
%       e.g., 's', 'i', '.3f', ...
%           Special column formats:
%           - 'rownr', displays current row number
%           - 'time', displays time since header was printed as HH:mm:ss
%           - 'sum{%.3e & %.3e}', mean+std summary, '&' represents '+-' 
%    varargin - name-value pairs
%        <'ColumnWidths',colWidths> - numeric, column width
%        <'SaveContent',saveContent> - logical, whether content should be saved
%        <'FileName',filename> - char, filename
%        <'fid',fid> - numeric, file identifier
%
% Outputs:
%    table - CORAtable object
%
% Example: 
%    % init table
%    hvalues = {'Heading 1','Heading 2'};
%    formats = {'s','s'};
%    table = CORAtable('single',hvalues,formats);
%    
%    % print table
%    table.printHeader()
%    table.printContentRow({'v11','v12'});
%    table.printContentRow({'v21','v22'});
%    table.printContentRow({'v31','v32'});
%    table.printFooter();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ASCIItable, CORAtable/printTable

% Authors:       Tobias Ladner
% Written:       19-September-2024
% Last update:   05-March-2025 (TL, printTable)
%                15-May-2025 (TL, added 'rownr' and 'time' format)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


properties
    % main properties
    design
end

methods
    function table = CORAtable(design,hvalues,formats,varargin)
        % 1. parse input
        [design,hvalues,formats,tableoptions] = aux_parseInput(design,hvalues,formats);

        % 2. call constructor
        table = table@ASCIItable(hvalues,formats,tableoptions{:},varargin{:});

        % 3. save properties
        table.design = design;
    end
end

methods (Static)
    function printTable(matlabtable,varargin)
        % prints a Matlab table in the specified CORA design
        % Syntax:
        %    CORAtable.printTable(matlabtable)
        %    CORAtable.printTable(matlabtable,design,varargin)
        %    CORAtable.printTable(matlabtable,design,formats,varargin)
        % 
        % Inputs: 
        %    matlabtable - table
        %    all other parameters as for the CORAtable constructor
        %
        % See also: CORAtable
        
        % parse input
        inputArgsCheck({{matlabtable,'att','table'}});
        hvalues = matlabtable.Properties.VariableNames;

        % parse table to have atomic values
        matlabtable = aux_parseAtomicValues(matlabtable);

        % detect design and format
        design = setDefaultValues({'single'},varargin);
        [formats,varargin] = aux_detectFormats(matlabtable,varargin(2:end));

        % get column widths
        colWidths = aux_determineColumnWidths(matlabtable,hvalues,formats);

        % init table
        coratable = CORAtable(design,hvalues,formats,'ColumnWidths',colWidths,varargin{:});

        % print
       coratable.printHeader()
       for i=1:height(matlabtable)
           cvalues = arrayfun(@(j) matlabtable{i,j}, 1:numel(hvalues),'UniformOutput',false);    
           coratable.printContentRow(cvalues);
       end
       coratable.printFooter();

    end
end

methods (Access=protected)

    function [hvalues,colWidths] = formatHeaderAndComputeColWidths(table,varargin)
        t = table;
        switch t.design
            case 'latex'
                % requires rewriting of header values

                hvalues = cell(size(t.hvalues));
                % Wrap header values with \text{...} to avoid math-mode;
                % for summary columns add \multicolumn{...}.
                for i=1:length(t.hvalues)
                    % read out hvalue
                    hvalue = t.hvalues{i};
                    % escape '#'
                    hvalue = strrep(hvalue,'#','\#');

                    format = t.formats{i};
                    if startsWith(format,'sum')
                        hvalues{i} = sprintf(...
                            '\\multicolumn{2}{c}{\\text{%s}}',hvalue);
                    else
                        hvalues{i} = sprintf('\\text{%s}',hvalue);
                    end
                end

                % Adapt column width for nice printing; \text{...} and
                % \multicolumn{...} change the width of the headings.
                colWidths = max(t.colWidths,cellfun(@strlength,hvalues));

                % call super function for proper formatting
                [hvalues,colWidths] = formatHeaderAndComputeColWidths@ASCIItable(t,hvalues,colWidths);
            otherwise
                % use default function defined in super class
                [hvalues,colWidths] = formatHeaderAndComputeColWidths@ASCIItable(t);
        end 
    end

    function format = getFormatString(table,format,colWidth) 
        if startsWith(format,'sum')
            % Replace separation indicator '&' with '±' (if design is not 
            % 'latex'; otherwise the latex table automatically inserts 
            % '\pm').
            if ~strcmp(table.design,'latex')
                format = strrep(format,'&','±');
            end
        end
        % Call super-class.
        format = getFormatString@ASCIItable(table,format,colWidth);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [design,hvalues,formats,tableoptions] = aux_parseInput(design,hvalues,formats)

    % check input args
    inputArgsCheck({ ...
        {design,'str',{'single', 'double', 'modern', 'minimalistic','ascii','latex','html','markdown','csv'}};    ...
        {hvalues,'att','cell'};    ...
        {formats,'att','cell'};    ...
    })

    % load design present
    switch design
        case 'single'
            tableoptions = aux_getSingleTableOptions();
        case 'double'
            tableoptions = aux_getDoubleTableOptions();
        case 'modern'
            tableoptions = aux_getModernTableOptions();
        case 'minimalistic'
            tableoptions = aux_getMinimalisticTableOptions();
        case 'ascii'
            tableoptions = aux_getAsciiTableOptions();
        case 'latex'
            tableoptions = aux_getLatexTableOptions(formats);
        case 'html'
            tableoptions = aux_getHTMLTableOptions();
        case 'markdown'
            tableoptions = aux_getMarkdownTableOptions();
        case 'csv'
            tableoptions = aux_getCSVTableOptions();
        otherwise
            % should have been caught at inputArgsCheck
            throw(CORAerror('CORA:wrongValue','third','Unknown design.'))
    end
end

function tableoptions = aux_getSingleTableOptions()
    % design: 'single' 
    tableoptions = {
        'tbhline','-','tbhcorner','|','hbvline','|','hsep','|','csep','|'
    };
end

function tableoptions = aux_getDoubleTableOptions()
    % design: 'double' 
    tableoptions = {
        'tbhline','=','tbhcorner','‖','hbvline','‖','hsep','|'
    };
end

function tableoptions = aux_getModernTableOptions()
    % design: 'modern' 
    tableoptions = {
        'tbhline','_','tbhcorner','','mbhline','-','mbhcorner','','bbhline','_','bbhcorner','', ...
        'hbvline','','hsep','|'
    };
end

function tableoptions = aux_getMinimalisticTableOptions()
    % design: 'minimalistic' 
    tableoptions = {
        'tbhline','-','tbhcorner','','hbvline','','hsep',' '
    };
end

function tableoptions = aux_getAsciiTableOptions()
    % design: 'ascii' 
    tableoptions = {
        'tbhline','-','tbhcorner','+','tbhsep','+','mbhsep','+','bbhsep','+','hbvline','|','hsep','|'
    };
end

function tableoptions = aux_getLatexTableOptions(formats)
    % design: 'latex' 

    % pre and post table text
    tpre = sprintf([ ...
        '%%%% Move to preamble.\n' ...
        '%%%% \\\\usepackage{booktabs} %%%% requires ''booktabs'' package\n' ...
        '%%%% \\\\usepackage{amstext} %%%% for \\\\text macro in headers\n' ...
        '%%%% \\\\usepackage{array} %%%% for \\\\newcolumntype macro\n' ...
        '%%%% \\\\newcolumntype{L}{>{$}l<{$}} %%%% math-mode version of "l" column type\n'...
        '%%%% \\\\newcolumntype{R}{>{$}r<{$}} %%%% math-mode version of "r" column type\n\n'...
        '\\\\begin{table}\n' ...
        '  \\\\centering\n' ...
        '  \\\\caption{My caption}\n' ...
        '  \\\\label{tab:my-label}\n' ...
        '  \\\\begin{tabular}{%s}\n' ...
    ], strjoin(cellfun(@aux_format2latex, formats,'UniformOutput',false),' '));
    tpost = [ ...
        '  \\end{tabular}\n' ...
        '\\end{table}\n' ...
    ];

    % build table options
    tableoptions = { ...
        'tbhline',' ','tbhcorner','','hbvline','','hsep','&', ...
        'tlpre','    \toprule','mlpre','    \midrule','blpre','    \bottomrule','hlpost','\\' ...
        'hlpre', '   ', 'tpre', tpre, 'tpost', tpost, ...
    };
end

function tableoptions = aux_getHTMLTableOptions(formats)
    % design: 'html' 
    tpre = '<table>\n';
    tpost = '</table>\n';

    tableoptions = { ...
        'tbhline',' ','tbhcorner','','hbvline','', ...
        'hlpre', '', 'tpre', tpre, 'tpost', tpost, ...
        'tlpre','  <thead>','mlpre','  </thead><tbody>','blpre','  </tbody>', ...
        'hlpre','    <tr><th>','hsep','</th><th>','hlpost','</th></tr>', ...
        'clpre','    <tr><td>','csep','</td><td>','clpost','</td></tr>', ...
    };
end

function tableoptions = aux_getMarkdownTableOptions()
    % design: 'markdown' 
    tableoptions = { ...
        'tbhline','','tbhcorner','','hbvline','|', ...
        'mbhline','-','mbhcorner','|','mbhsep','|'
    };
end

function tableoptions = aux_getCSVTableOptions()
    % design: 'csv' 
    tableoptions = {
        'tbhline','','tbhcorner','','hbvline','','hsep',';'
    };
end

function latexformat = aux_format2latex(format)
    if startsWith(format,'sum')
        % Insert \pm in between.
        latexformat = 'R@{$\\pm$}L';
    elseif contains(format,'s')
        % align strings on the left
        latexformat = 'l'; 
    else
        % align numbers etc. on the right; format in math-mode
        latexformat = 'R';
    end
end

function [formats,givenvalues] = aux_detectFormats(matlabtable,givenvalues)

    % check if formats are given
    if numel(givenvalues) > 1 && iscell(givenvalues{1})
        % assume first entry is formats
        formats = givenvalues{1};
        givenvalues = givenvalues(2:end);
        return
    end

    % try to guess the format from given table ---

    % extract variable types
    types = aux_getColumnTypes(matlabtable);
    
    % determine formats
    formats = cell(1,numel(types));
    for i=1:numel(types)
        type = types(i);
        if ismember(type,{'double','logical'})
            % numbers are printed with %d option
            formats{i} = 'd';
        else
            % default as string
            formats{i} = 's';
        end
    end
end

function colWidths = aux_determineColumnWidths(matlabtable,hvalues,formats)
    % init
    colWidths = cellfun(@(hvalue) strlength(hvalue),hvalues,'UniformOutput',true);

    % determine column widths
    for j=1:numel(formats)
        format = formats{j};
        % add '%' if necessary
        if ~contains(format,'%')
            format = sprintf('%%%s',format);
        end

        % determine widths of header and all entries
        widths = [ ...
            colWidths(j), ...
            arrayfun(@(i) aux_determineLengthEntry(matlabtable{i,j},format), 1:height(matlabtable)) ...
        ];

        % save max width
        colWidths(j) = max(widths);
    end
end

function len = aux_determineLengthEntry(entry,format)
    % format entry
    formattedEntry = sprintf(format,entry);
    % determine length (ignoring hyperlink text)
    [~,len] = centerString(formattedEntry,0);
end

function matlabtable = aux_parseAtomicValues(matlabtable)
    % parses the matlab table to only have atomic values
    
    % find cell columns
    types = aux_getColumnTypes(matlabtable);
    idx = find(strcmp('cell',types));    
    for i=idx
        % convert each column

        % read out column
        column = matlabtable(:,i);
        columnname = column.Properties.VariableNames{1};

        formatfun = @(value) sprintf("%s",value{1});
        if height(column) > 1
            value = column{1,1};
            if ismember(class(value{1}),{'double','logical'})
                formatfun = @(value) value{1};
            end
        end

        % remove outer cell
        column = rowfun(formatfun, column);
        column.Properties.VariableNames{1} = columnname;

        % store back into table
        matlabtable = [matlabtable(:,1:(i-1)),column,matlabtable(:,(i+1):end)];

    end

end

function types = aux_getColumnTypes(matlabtable)
    if height(matlabtable) >= 1
        if isMATLABReleaseOlderThan('R2024b')
            % determine class based on first row
            types = arrayfun(@(i) class(matlabtable{1,i}),1:numel(matlabtable.Properties.VariableNames),'UniformOutput',false);
        else
            % select variable types stored in table
            types = matlabtable.Properties.VariableTypes;
        end
    else
        % just use string; no entries anyway
        types = arrayfun(@(i) 'string',1:numel(matlabtable.Properties.VariableNames),'UniformOutput',false);
    end
end

% ------------------------------ END OF CODE ------------------------------
