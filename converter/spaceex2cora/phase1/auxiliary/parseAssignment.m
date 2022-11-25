function [varnames,exprs,warnings] = parseAssignment(str)
%	INTPUTS:
%   str (string): series of assignments seperated by "&"
%               assignment syntax: <variable>' = <expr>
%               operators = == := supported
%   OUTPUTS:
%   varnames (string array): names of assigned variables
%   exprs (symbolic array): assigned expressions
%       represent assignments "varnames(i)' = exprs(i)"
%   warnings (struct) document failed operations in warnings.messages
 warnings = struct([]);
 warn_ct = 0;
 
 % Unfortunately, the symbolic Toolbox interprets variables named "i","j",
 % "I" or "J" as the imaginary number.
 % We perform a transformation on all variable names to avoid this.
 str = replaceImagVarnames(str);
 
 % if char array was passed, transform to string
 str = string(str);
 
 % split into single equations at '&' signs
 equations = strsplit(str,"&");
 
 % preallocate string arrays for left & right sides
 varnames = strings(length(equations),1);
 exprStrings = strings(length(equations),1);
 numExpr = 0;
 
 % split equations into sides
 for i = 1:length(equations)
    % split by assignment operators =, ==, or :=
    sides = strsplit(equations(i),"(=|:)?=","DelimiterType","RegularExpression");
    if length(sides) == 2
        
        % trim variable name, then remove uptick by going back to char array
        temp = strtrim(sides(1));
        temp = char(temp);
        if(temp(end) == '''')
            % trim uptick
            temp = temp(1:end-1);
        else
            warn_ct = warn_ct+1;
            warnings(warn_ct).message = sprintf(...
                "Format Warning: assigned variable %s has no uptick",...
                temp);
        end
        
        % update counter & store processed strings
        numExpr = numExpr+1;
        varnames(numExpr) = string(temp);
        exprStrings(numExpr) = sides(2);
    else
        warn_ct = warn_ct+1;
        warnings(warn_ct).message = sprintf(...
            "Format Error: %i operators found in assignment ""%s"", skipping",...
            length(sides)-1,equations(i));
    end
 end
 
 % shorten arrays, in case of invalid equations
 varnames = varnames(1:numExpr,1);
 % convert expressions to symbolic
 exprs = str2symbolic(exprStrings(1:numExpr,1));
end