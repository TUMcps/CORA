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
%    formats - cell array containing cell formats; for summary columns
%       'sum{%.3e & %.3e}', where '&' indicates the position of '+-' 
%       (or \pm in latex design).
%    varargin - name-value pairs
%        <'ColumnWidths',colWidths> - numeric, column width
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
% See also: ASCIItable

% Authors:       Tobias Ladner
% Written:       19-September-2024
% Last update:   ---
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

methods (Access=protected)

    function [hvalues,colWidths] = formatHeaderAndComputeColWidths(table,varargin)
        t = table;
        switch t.design
            case 'latex'
                hvalues = cell(size(t.hvalues));
                % Wrap header values with \text{...} to avoid math-mode;
                % for summary columns add \multicolumn{...}.
                for i=1:length(t.hvalues)
                    format = t.formats{i};
                    if startsWith(format,'sum')
                        hvalues{i} = sprintf(...
                            '\\multicolumn{2}{c}{\\text{%s}}',t.hvalues{i});
                    else
                        hvalues{i} = sprintf('\\text{%s}',t.hvalues{i});
                    end
                end

                % Adapt column width for nice printing; \text{...} and
                % \multicolumn{...} change the width of the headings.
                colWidths = max(t.colWidths,cellfun(@strlength,hvalues));

                % call super function for proper formatting
                [hvalues,colWidths] = ...
                    t.formatHeaderAndComputeColWidths@ASCIItable(hvalues,colWidths);
            otherwise
                [hvalues,colWidths] = ...
                    t.formatHeaderAndComputeColWidths@ASCIItable();
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
        format = table.getFormatString@ASCIItable(format,colWidth);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [design,hvalues,formats,tableoptions] = aux_parseInput(design,hvalues,formats)

    % check input args
    inputArgsCheck({ ...
        {design,'str',{'single', 'double', 'modern', 'minimalistic','ascii','latex','html','markdown'}};    ...
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
        otherwise
            % should have been caught at inputArgsCheck
            throw(CORAerror('CORA:wrongValue','third','Unknown design.'))
    end
end

function tableoptions = aux_getSingleTableOptions()
    tableoptions = {
        'tbhline','-','tbhcorner','|','hbvline','|','hsep','|','csep','|'
    };
end

function tableoptions = aux_getDoubleTableOptions()
    tableoptions = {
        'tbhline','=','tbhcorner','‖','hbvline','‖','hsep','|'
    };
end

function tableoptions = aux_getModernTableOptions()
    tableoptions = {
        'tbhline','_','tbhcorner','','mbhline','-','mbhcorner','','bbhline','_','bbhcorner','', ...
        'hbvline','','hsep','|'
    };
end

function tableoptions = aux_getMinimalisticTableOptions()
    tableoptions = {
        'tbhline','-','tbhcorner','','hbvline','','hsep',' '
    };
end

function tableoptions = aux_getAsciiTableOptions()
    tableoptions = {
        'tbhline','-','tbhcorner','+','tbhsep','+','mbhsep','+','bbhsep','+','hbvline','|','hsep','|'
    };
end

function tableoptions = aux_getLatexTableOptions(formats)
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

    tableoptions = { ...
        'tbhline',' ','tbhcorner','','hbvline','','hsep','&', ...
        'tlpre','    \toprule','mlpre','    \midrule','blpre','    \bottomrule','hlpost','\\' ...
        'hlpre', '   ', 'tpre', tpre, 'tpost', tpost, ...
    };
end

function tableoptions = aux_getHTMLTableOptions(formats)
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
    tableoptions = { ...
        'tbhline','','tbhcorner','','hbvline','|', ...
        'mbhline','-','mbhcorner','|','mbhsep','|'
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

% ------------------------------ END OF CODE ------------------------------
