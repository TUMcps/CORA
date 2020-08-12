function safetext = replaceImagVarnames(text)
% replaceImagVarnames - given a string mathematical expression, replace all
%    variable names, that could interpreted as the imaginary number
%    by the symbolic toolkit.
%    Replace (i,j,I,J) by safe names (ii,jj,II,JJ).
%    Apply the following name transformation to conserve name uniqueness:
%       transform(str) = c^(n+1)  IF  str = c^n & c ∈ {i,j,I,J}
%       transform(str) = str        OTHERWISE
%
% Syntax:  
%    safetext = replaceImagVarnames(text)
%
% Inputs:
%    text - expression in text (string nxm OR char 1xn)
%
% Outputs:
%    safetext - expression with same logic as text, safe for sym parsing
%
% Example: replaceImagVarnames("I + 5*(time - jj)") -> "II + 5*(time - jjj)"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Johann Schöpfer
% Written: 09-April-2018 
% Last update: 09-April-2018 
% Last revision: ---

%------------- BEGIN CODE --------------

% FOR MULTIPLE STRINGS: self-call this function for each individual string
if ~ischar(text) && ~isStringScalar(text)
    safetext = strings(size(text));
    for i = 1:numel(text)
        safetext(i) = replaceImagVarnames(text(i));
    end
    
    return;
end

% FOR SCALAR STRINGS OR CHAR VECTORS: perform char manipulations

% This regular expression detecs variable names in arithmetic expressions,
% by finding terms surrounded by non-variable characters.
% Then find variables that match any of these patterns: i+ j+ I+ J+
regex = '(?<!\w)(i+|j+|I+|J+)(?!\w)';

% run regexp on char vector for easier manipulation later
chars = char(text);

% evaluate regex, store matches and their start location
[matches,idxs] = regexp(chars,regex,'match','start');

% reassemble the text in safechars
safechars = [];
% apply the transformation by duplicating the first letter of each match
% index the first char in chars, that has NOT been copied to the new array
lastIdx = 1;
% itearate over matches
for i = 1:numel(matches)
    matchStart = idxs(i);
    firstLetter = matches{i}(1);

    % copy all chars leading up to this match
    safechars = [safechars chars(lastIdx:matchStart-1)];
    
    % add the additional letter for the current match
    safechars = [safechars firstLetter];
    
    % update the pointer for not yet copied chars
    lastIdx = matchStart;
end
% finally, copy any leftover text
safechars = [safechars chars(lastIdx:end)];

% convert the char vector back to string, & return the result
safetext = string(safechars);

end

%------------- END OF CODE -------------