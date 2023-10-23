function [inequalities,equalities,exprLeft,exprRight,warnings] = parseCondition(str)
% parseCondition - parses a string of equalities and inequalities (aka the
%    conditions) which define an invariant or a guard set and returns
%    symbolic expressions
%
% Syntax:
%    [inequalities,equalities,exprLeft,exprRight,warnings] = parseCondition(str)
%
% Inputs:
%    str (string) - series of [in-]equalities separated by '&'
%                   operators < <= > >= = == := supported
%                   (</<= & >/>= are treated identically)
%                   chain [in-]equalities (i.e. 'w < x = y < z') allowed
%
% Outputs:
%    inequalities (symbolic array) - each sym represents inequality "sym <= 0"
%    equalities (symbolic array) - each sym represents equality "sym == 0"
%    exprLeft - ??? (Note: seems to only be relevant if sx model is
%                    conversed to strictly flat hybrid automaton)
%    exprRight - ???
%    warnings (struct) - document failed operations in warnings.messages
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

warnings = struct([]);
warn_ct = 0;

% Unfortunately, the symbolic Toolbox interprets variables named "i","j",
% "I" or "J" as the imaginary number.
% We perform a transformation on all variable names to avoid this.
str = replaceImagVarnames(str);

% if char array was passed, transform to string
str = string(str);

% split into single equations (delimiter '&' is mandatory in SpaceEx)
equations = strsplit(str,'&');

% split equations by side
% reassign left/right side s.t. "left <= right" for inequalities
% doesn't matter for equalities

% preallocate arrays with estimate size
sz = [length(equations),1];
IneqExprsLeft = strings(sz);
IneqExprsRight = strings(sz);
EqExprsLeft = strings(sz);
EqExprsRight = strings(sz);
% count found conditions
numIneq = 0;
numEq = 0;

for i = 1:length(equations)

    [sides,operators] = strsplit(equations(i),'((=|:)?=)|((<|>)=?)',...
        'DelimiterType','RegularExpression');
    
    for j = 1:(length(sides)-1)
        if length(strfind(operators{j},'<')) == 1
            % add "less than" inequality
            numIneq = numIneq + 1;
            IneqExprsLeft(numIneq,1) = sides(j);
            IneqExprsRight(numIneq,1) = sides(j+1);
        elseif length(strfind(operators{j},'>')) == 1
            % add "greater than" inequality
            % (switch sides resp. "less than")
            numIneq = numIneq + 1;
            IneqExprsLeft(numIneq,1) = sides(j+1);
            IneqExprsRight(numIneq,1) = sides(j);
        else
            % add equality
            numEq = numEq + 1;
            EqExprsLeft(numEq,1) = sides(j);
            EqExprsRight(numEq,1) = sides(j+1);
        end
    end

end

% bring expressions into form "sym <= 0" or "sym == 0", respectively, by
% converting sides to symbolic, then subtracting right from left sides
inequalities = str2symbolic(IneqExprsLeft(1:numIneq,1))...
                      - str2symbolic(IneqExprsRight(1:numIneq,1));
equalities = str2symbolic(EqExprsLeft(1:numEq,1))...
                    - str2symbolic(EqExprsRight(1:numEq,1));

% added: (also to return list)
exprLeft = str2symbolic(EqExprsLeft(1:numEq,1));
exprRight = str2symbolic(EqExprsRight(1:numEq,1));

% ------------------------------ END OF CODE ------------------------------
