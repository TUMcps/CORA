function exportScriptToHTML(filepath)
% exportScriptToHTML - exports the given script to html
%
% Syntax:
%    exportScriptToHTML(filepath)
%
% Inputs:
%    filepath - char, file path
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       18-July-2023
% Last update:   30-October-2024 (TL, generalized)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,1);

% turn off warnings
w = warning;
warning off

% get file extension
[~,~,ext] = fileparts(filepath);

% convert script to html
htmlpath = convertStringsToChars(strrep(filepath, ext, ".html"));
matlab.internal.liveeditor.openAndConvert(filepath, htmlpath);

% remove unnecessary padding ---

% read file
filetext = fileread(htmlpath);

% adjust css
filetext = strrep(filetext,'.rtcContent { padding: 30px;', '.rtcContent { padding: 0px;');
filetext = strrep(filetext,'.CodeBlock { ', '.CodeBlock { width: 100%; ');

% replace all % with %% for fprintf
filetext = strrep(filetext,'%', '%%');

% write to file
fid = fopen(htmlpath,'w');
fprintf(fid,filetext);
fclose(fid);

% restore warning
warning(w)

end

% ------------------------------ END OF CODE ------------------------------
