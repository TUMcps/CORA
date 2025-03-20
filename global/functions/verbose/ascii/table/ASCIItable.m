classdef (Abstract) ASCIItable < handle
% ASCIItable - creates a nice ASCII table
%
% Syntax:
%    table = ASCIItable(hvalues,formats,tableoptions)
%
% Inputs:
%    hvalues - cell array containing headings
%    formats - cell array containing cell formats
%    tableoptions - name-value pair
%        <'ColumnWidths',colWidths> - numeric, column widths
%        <'SaveContent',saveContent> - logical, whether content should be saved
%        ------------------------------------------------------------------
%        <'tbhline',tbhline> - char, top boundary horizontal line
%        <'tbhcorner',tbhcorner> - char, top heading boundary corner
%        <'tbhsep',tbhsep> - char, top heading boundary separator
%        <'mbhline',tbhline> - char, mid boundary horizontal line
%        <'mbhcorner',tbhcorner> - char, mid heading boundary corner
%        <'mbhsep',mbhsep> - char, mid heading boundary separator
%        <'bbhline',tbhline> - char, bottom boundary horizontal line
%        <'bbhcorner',tbhcorner> - char, bottom heading boundary corner
%        <'bbhsep',bbhsep> - char, bottom heading boundary separator
%        <'hbvline',hbvline> - char, heading boundary vertical line
%        <'hsep',hsep> - char, heading separator
%        <'cbvline',cbvline> - content boundary vertical line
%        <'csep',csep> - char, content separator
%        <'tpre',tbegin> - char, pre table text
%        <'tpost',tbegin> - char, post table text
%        <'tlpre',tbegin> - char, pre top line text
%        <'mlpre',tbegin> - char, pre mid line text
%        <'blpre',tbegin> - char, pre bottom line text
%        <'hlpre',tbegin> - char, pre header line text
%        <'clpre',tbegin> - char, pre content line text
%        <'tlpost',tbegin> - char, post top line text
%        <'mlpost',tbegin> - char, post mid line text
%        <'blpost',tbegin> - char, post bottom line text
%        <'hlpost',tbegin> - char, post header line text
%        <'clpost',tbegin> - char, post content line text
%        ------------------------------------------------------------------
%
% Outputs:
%    table - ASCIItable object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAtable

% Authors:       Lukas Koller, Tobias Ladner
% Written:       19-September-2024
% Last update:   05-March-2025 (TL, option to save content)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


properties
    % main properties
    hvalues
    formats
    colWidths
    saveContent
    content = {};
end

properties (Hidden)
    % boundary chars
    tbhline,tbhcorner,tbhsep
    mbhline,mbhcorner,mbhsep
    bbhline,bbhcorner,bbhsep
    % heading chars
    hbvline,hsep
    % content chars
    cbvline,csep
    % pre/post
    tpre,tpost
    tlpre,tlpost
    hlpre,hlpost
    mlpre,mlpost
    blpre,blpost
    clpre,clpost
end

methods
    function table = ASCIItable(hvalues,formats,varargin)
        % 1. parse input
        [hvalues,formats,colWidths,saveContent, ...
            tbhline,tbhcorner,tbhsep, ...
            mbhline,mbhcorner,mbhsep, ...
            bbhline,bbhcorner,bbhsep, ...
            hbvline,hsep, ...
            cbvline,csep, ...
            tpre,tpost, ...
            tlpre,tlpost, ...
            hlpre,hlpost, ...
            mlpre,mlpost, ...
            blpre,blpost, ...
            clpre,clpost] = aux_parseInput(hvalues,formats,varargin);

        % 2. save properties
        table.hvalues = hvalues;
        table.formats = formats;
        table.colWidths = colWidths;
        table.saveContent = saveContent;
        % boundary chars
        table.tbhline = tbhline;
        table.tbhcorner = tbhcorner;
        table.tbhsep = tbhsep;
        table.mbhline = mbhline;
        table.mbhcorner = mbhcorner;
        table.mbhsep = mbhsep;
        table.bbhline = bbhline;
        table.bbhcorner = bbhcorner;
        table.bbhsep = bbhsep;
        % heading chars
        table.hbvline = hbvline;
        table.hsep = hsep;
        % heading chars
        table.cbvline = cbvline;
        table.csep = csep;
        % pre/post
        table.tpre = tpre;
        table.tpost = tpost;
        table.tlpre = tlpre;
        table.tlpost = tlpost;
        table.hlpre = hlpre;
        table.hlpost = hlpost;
        table.mlpre = mlpre;
        table.mlpost = mlpost;
        table.blpre = blpre;
        table.blpost = blpost;
        table.clpre = clpre;
        table.clpost = clpost;
    end

    function printHeader(table)
        % print heading block
        t = table;
        
        % Determine column widths and format header values; e.g. wrap in 
        % \text{...} for latex design to avoid math-mode in CORAtable.
        [hvalues,colWidths] = formatHeaderAndComputeColWidths(t);

        % build heading row
        hrow = [...
          t.hlpre t.hbvline ' ' ...
          strjoin(hvalues, [' ' t.hsep ' ']) ...
          ' ' t.hbvline t.hlpost ...
        ];

        % print
        fprintf(t.tpre)
        table.printTopBoundaryRow();
        disp(hrow);
        table.printMidBoundaryRow();
    end
 
    function printContentRow(table,cvalues)
        % print content row

        % Format column values.
        fcolumns = formatColumnValues(table,cvalues);

        % build content row
        t = table;
        crow = [ ...
            ... % start row with
            t.clpre t.cbvline ' ' ...
            ... % join all cells
            strjoin(fcolumns,[' ' t.csep ' ']) ...
            ... % and finish with
            ' ' t.cbvline t.clpost ...
        ];

        % print
        disp(crow);

        % save
        if table.saveContent
            table.content(end+1,:) = cvalues;
        end
    end
    
    function printFooter(table)
        % print bottom of table
        table.printBottomBoundaryRow();
        fprintf(table.tpost)
        disp(' ')
    end

    % helper ---

    function printTopBoundaryRow(table)
        % print top boundary row "\toprule"
        disp(buildBoundaryRow(table,table.tbhline,table.tbhcorner, ...
            table.tbhsep,table.tlpre,table.tlpost));
    end

    function printMidBoundaryRow(table)
        % print top boundary row "\midrule"
        disp(buildBoundaryRow(table,table.mbhline,table.mbhcorner, ...
            table.mbhsep,table.mlpre,table.mlpost));
    end

    function printBottomBoundaryRow(table)
        % print top boundary row "\bottomrule"
        disp(buildBoundaryRow(table,table.bbhline,table.bbhcorner, ...
            table.bbhsep,table.blpre,table.blpost));
    end

    function reprint(table,varargin)
        % reprints the saved table as CORAtable with the specified design
        
        % parse input
        design = setDefaultValues('single',varargin);

        % show warning if no content was saved
        if isempty(table.content)
            CORAwarning("CORA:global",'Table has no saved content. Make sure to set ''SaveContent'' to true when initializing the table.')
        end

        % reprint table
        coratable = CORAtable(design,table.hvalues,table.formats,"ColumnWidths",table.colWidths);
        coratable.printHeader();
        for i=1:size(table.content,1)
            coratable.printContentRow(table.content(i,:))
        end
        coratable.printFooter();
    end
end

methods (Access=protected)

    function [hvalues,colWidths] = formatHeaderAndComputeColWidths(table,hvalues,colWidths)
        % format header text and determine column width 

        % parse input
        if nargin < 3
            hvalues = table.hvalues;
            colWidths = table.colWidths;
        end

        % pre/append headings to match column widths
        hvalues = arrayfun( ...
            @(i) centerString(hvalues{i},colWidths(i)), ...
            1:numel(hvalues),'UniformOutput',false);
    end

    function fcolumns = formatColumnValues(table,cvalues)
        t = table;
        
        % Compute column widths, which depend on the design 
        % (see CORAtable); thus are computed dynamically.
        [~,colWidths] = formatHeaderAndComputeColWidths(table);
        
        % Store summarized values.
        fcolumns = {};
        % Compute statistics for each summary column.
        for i=1:length(cvalues)
            % Obtain format string.
            format = getFormatString(t,t.formats{i},colWidths(i));
            % fill in values
            if startsWith(t.formats{i},'sum')
                formattedValues = sprintf(format,mean(cvalues{i}),std(cvalues{i}));
            else
                formattedValues = sprintf(format,cvalues{i});
            end
            % center values
            % (alignment is usually already handled in getFormatString,
            % but we re-center it here if the column-width changed in between)
            fcolumns{end+1} = centerString(formattedValues,colWidths(i));
        end
    end

    function format = getFormatString(table,format,colWidth)
        % There is an auxiliary function to avoid duplication of code; this
        % function is overwritten in CORAtable to exchange the separation
        % indicator for summary-cells.
        format = aux_getFormatString(format,colWidth);
    end

    function hbrow = buildBoundaryRow(table,bhline,bhcorner,bhsep,lpre,lpost)
        t = table;
        
        % Compute column widths, which depend on the design 
        % (see CORAtable); thus are computed dynamically.
        [~,colWidths] = formatHeaderAndComputeColWidths(t);

        % build heading boundary row
        hbrow = [ ...
            ... % start row with
            lpre bhcorner bhline ...
            ... % join all cells
            strjoin( ...
                ... % build boundary with correct cell length
                arrayfun(@(width) repmat(bhline,[1,width]),colWidths, ...
                    'UniformOutput',false), ...
                ... % separate each cell with
                [bhline bhsep bhline]) ...
            ... % and finish with
            bhline bhcorner lpost...
        ];
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [hvalues,formats,colWidths,saveContent, ...
            tbhline,tbhcorner,tbhsep, ...
            mbhline,mbhcorner,mbhsep, ...
            bbhline,bbhcorner,bbhsep, ...
            hbvline,hsep, ...
            cbvline,csep, ...
            tpre,tpost, ...
            tlpre,tlpost, ...
            hlpre,hlpost, ...
            mlpre,mlpost, ...
            blpre,blpost, ...
            clpre,clpost] = aux_parseInput(hvalues,formats,NVpairs)

    % parse input
    inputArgsCheck({ ...
        {hvalues,'att','cell','isrow'}; ...
        {formats,'att','cell','isrow'}; ...
    })

    % read name-value pairs
    [NVpairs,colWidths] = readNameValuePair(NVpairs,'ColumnWidths','isnumeric',0);
    [NVpairs,saveContent] = readNameValuePair(NVpairs,'SaveContent','islogical',false);
    % boundary chars
    [NVpairs,tbhline] = readNameValuePair(NVpairs,'tbhline','ischar','-');
    [NVpairs,tbhcorner] = readNameValuePair(NVpairs,'tbhcorner','ischar','+');
    [NVpairs,tbhsep] = readNameValuePair(NVpairs,'tbhsep','ischar',tbhline);
    [NVpairs,mbhline] = readNameValuePair(NVpairs,'mbhline','ischar',tbhline);
    [NVpairs,mbhcorner] = readNameValuePair(NVpairs,'mbhcorner','ischar',tbhcorner);
    [NVpairs,mbhsep] = readNameValuePair(NVpairs,'mbhsep','ischar',mbhline);
    [NVpairs,bbhline] = readNameValuePair(NVpairs,'bbhline','ischar',tbhline);
    [NVpairs,bbhcorner] = readNameValuePair(NVpairs,'bbhcorner','ischar',tbhcorner);
    [NVpairs,bbhsep] = readNameValuePair(NVpairs,'bbhsep','ischar',bbhline);
    % heading chars
    [NVpairs,hbvline] = readNameValuePair(NVpairs,'hbvline','ischar','|');
    [NVpairs,hsep] = readNameValuePair(NVpairs,'hsep','ischar','|');
    % content chars
    [NVpairs,cbvline] = readNameValuePair(NVpairs,'cbvline','ischar',hbvline);
    [NVpairs,csep] = readNameValuePair(NVpairs,'csep','ischar',hsep);
    % pre/post
    [NVpairs,tpre] = readNameValuePair(NVpairs,'tpre','ischar','');
    [NVpairs,tpost] = readNameValuePair(NVpairs,'tpost','ischar','');
    [NVpairs,tlpre] = readNameValuePair(NVpairs,'tlpre','ischar','');
    [NVpairs,tlpost] = readNameValuePair(NVpairs,'tlpost','ischar','');
    [NVpairs,hlpre] = readNameValuePair(NVpairs,'hlpre','ischar',tlpost);
    [NVpairs,hlpost] = readNameValuePair(NVpairs,'hlpost','ischar','');
    [NVpairs,mlpre] = readNameValuePair(NVpairs,'mlpre','ischar','');
    [NVpairs,mlpost] = readNameValuePair(NVpairs,'mlpost','ischar',tlpost);
    [NVpairs,blpre] = readNameValuePair(NVpairs,'blpre','ischar','');
    [NVpairs,blpost] = readNameValuePair(NVpairs,'blpost','ischar',tlpost);
    [NVpairs,clpre] = readNameValuePair(NVpairs,'clpre','ischar',hlpre);
    [NVpairs,clpost] = readNameValuePair(NVpairs,'clpost','ischar',hlpost);

    % check name-value pairs

    % correctly set column widths
    if isscalar(colWidths)
        colWidths = ones(size(hvalues)) * colWidths;
    end

    % Compute width of formatted columns.
    formatWidth = zeros(size(colWidths));
    numColElems = cellfun(@(f) max(1,sum(f == '%')),formats);
    for i=1:length(formats)
        % Create enough format arguments.
        args = repmat({0},1,numColElems(i));
        % Extract format strings.
        format = aux_getFormatString(formats{i},0);
        % Compute length of a formatted string.
        formatWidth(i) = strlength(sprintf(format,args{:}));
    end

    % Make columns width at least as long as headings and format widths.
    colWidths = max([colWidths; cellfun(@strlength,hvalues); formatWidth]);

    if CHECKS_ENABLED
        % length of hvalues and formats must match
        if numel(hvalues) ~= numel(formats)
            throw(CORAerror('CORA:wrongValue','first/second','Length of headings and formats must match.'))
        end

        % check if no format has '%' (gets added in printContentRow to determine cell width)
        if any(cellfun(@(format) contains(format,'%') && ~startsWith(format,'sum'),formats,'UniformOutput',true))
            throw(CORAerror('CORA:wrongValue','second','Specify formats without leading ''%%''. This is required to set the cell width correctly.'))
        end

        % % length of hvalues and colWidths must match
        % if numel(hvalues) ~= numel(colWidths)
        %     throw(CORAerror('CORA:wrongValue','first/second','Length of headings and column widths must match.'))
        % end
    end
end

function format = aux_getFormatString(format,colWidth)
    if startsWith(format,'sum')
        % Remove enclosing 'sum{...}'.
        format = regexp(format,'(?<=sum{).+(?=})','match','once');
    else
        % Add '%' and ensure proper length of formatted string.
        if colWidth > 0
            if strcmp(format,'s')
                % strings are left aligned
                format = ['-' num2str(colWidth) format];
            else
                % everything else is right aligned
                format = [num2str(colWidth) format];
            end
        end
        format = ['%' format];
    end
end

% ------------------------------ END OF CODE ------------------------------
