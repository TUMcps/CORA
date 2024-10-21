function display(obj)
% display - Display the STL interval to the command window
%
% Syntax:
%    display(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    ---
%
% Example:
%    I = stlInterval(0,1,true,false);
%    display(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);
disp([inputname(1), ' =']);
fprintf(newline);

if ~isscalar(obj)
    sizeStr = strrep(mat2str(size(obj)), ' ', 'x');
    sizeStr = sizeStr(2:end-1); % remove brackets
    fprintf('  %s <a href="matlab:helpPopup %s">%s</a> array\n\n',sizeStr,class(obj),class(obj));
end
for i = 1:numel(obj)
    fprintf('  %s\n\n',obj(i).toStr());
end

% ------------------------------ END OF CODE ------------------------------
