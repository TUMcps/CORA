function str = printMatrixConverter(M,maxLineLength)
% print matrix M in compact matlab-readable format
% prints numbers with 32 digits of precision

% 2nd function: ensure result string has no exessively long lines
% split into multiple lines via ... if necessary

if(nargin<2)
    maxLineLength = 75;
end

% significant digits of printed numerals
precision = 32;

% sprintf format code for numerals with set precision(i.e. '%.10G')
numFormat = ['%.' int2str(precision) 'G'];

dims = size(M);

if ~all(dims)
    % retain dimensions of matrices w/o values
    str = sprintf("zeros(%s)",printMatrixConverter(dims));
elseif length(dims)>2
    % not supported yet
    str = "[]";
    warning("cannot print matrix of dimensions [%s], return []",num2str(dims));
else
    % begin array
    str = "[";
    
    % iterate over rows
    for i = 1:dims(1)
        if(i ~= 1)
            % add row separator
            str = str + ";";
        end
        
        % convert i-th row to comma-separated string
        row = sprintf([numFormat ','], M(i,:));
        
        % remove trailing comma & convert to string
        row = string(row(1:end-1));
        
        str = str + row;
    end
    % end matrix
    str = str + "]";
end

% insert linebreaks (with ...) if necessary
if strlength(str) > maxLineLength
   % convert to char array for easier splicing
   chars = char(str);
   str = "";
   
   % insert "..." and linebreak only after ',' or ';'
   isSeparator = (chars == ',') | (chars == ';');
   
   % go through chars_in, sectioning off maximal lines
   num_chars = length(chars);
   chars_pointer = 1;
   
   % while the remainder does not fit in 1 line, loop
   while (chars_pointer - 1) < (num_chars - maxLineLength)
       % find the last possible separator in the next line,
       % such that there is still space for the '...'
       max_idx = (chars_pointer + maxLineLength - 1) - 3;
       sep_idx = find(isSeparator(chars_pointer:max_idx),1,'last');
       % sep_idx is local index in searched space, make it global in chars
       sep_idx = sep_idx + chars_pointer - 1;
       
       
       % add the line to output
       str = str + chars(chars_pointer:sep_idx) + "..." + newline;
       
       % move up pointer
       chars_pointer = sep_idx + 1;
   end
   
   % add last line to output
   str = str + chars(chars_pointer:end);
end

end
