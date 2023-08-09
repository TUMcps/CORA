function exportToHTML()
% exportToHTML - exports the current open live script (*.mlx) to html
%
% Syntax:  
%    exportToHTML()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Tobias Ladner
% Written:      18-July-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% turn off warnings
warnOrg = warning;
warning off

% convert mlx to html
filePath = matlab.desktop.editor.getActiveFilename;
htmlPath = convertStringsToChars(strrep(filePath, ".mlx", ".html"));
matlab.internal.liveeditor.openAndConvert(filePath, htmlPath);

% remove "exportHTML" text
filetext = fileread(htmlPath);
pattern = '<span >exportToHTML(); </span><span style="color: rgb(0, 128, 19);">% for website</span></span>';
filetext = strrep(filetext,pattern,'');
pattern = '.rtcContent { padding: 30px; }';
filetext = strrep(filetext,pattern,'');
fid  = fopen(htmlPath,'w');
fprintf(fid,'%s',filetext);
fclose(fid);

% restore warning
warning(warnOrg)

%------------- END OF CODE --------------

