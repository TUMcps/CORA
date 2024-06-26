function names_out = applyRenames(names_in,keys,renames)
% applyRenames - applies a list of name changes to a string array
%
% Syntax:
%    names_out = applyRenames(names_in,keys,renames)
%
% Inputs:
%    names_in (string) - strings to be renamed
%    renames (2xN string) - string renamings to be applied
%                           change from names in 1st column to names in
%                           2nd column (created by findRenames.m)
%
% Outputs:
%    names_out - names_in with applied renames
%                (renames not to be applied recursively)
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

names_out = names_in;
% iterate over names_in and find the rename
for i = 1:numel(names_in)
    % find the 1st rename applying to names_in(i)
    rename_idx = find(keys == names_in(i),1);
    if ~isempty(rename_idx)
        % if a rename is found, write it to output
        names_out(i) = renames(rename_idx);
    end
end

% ------------------------------ END OF CODE ------------------------------
