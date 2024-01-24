function convertToTikz(figpath,varargin)

% open figure
openfig(figpath);
box on;
drawnow;

% make tikz image
% please ask for access of this repo: https://gitlab.lrz.de/cora/matlab2tikz
% then download it and add it to the Matlab path
tikzpath = strrep(figpath, '.fig', '.tikz');
matlab2tikz(tikzpath, 'relativePrecision', 5e-3, ...
    'floatFormatX', '%.4f', 'floatFormatY', '%.4f', ...
    'showInfo',false,varargin{:});
aux_postTikz(tikzpath);

% close figure;
% close;

end


% Auxiliary functions -----------------------------------------------------

function aux_postTikz(tikzpath)

fid  = fopen(tikzpath,'r');
f=fread(fid,'*char')';
fclose(fid);

f = regexprep(f,'width=[0-9]*.[0-9]*in','width=4cm');
f = regexprep(f,'height=[0-9]*.[0-9]*in','height=4cm');
f = regexprep(f,'at=\{\([0-9]*.[0-9]*in,[0-9]*.[0-9]*in\)\},','at={(0in,0in)},');
f = strrep(f,"\begin{tikzpicture}",compose("\\begin{tikzpicture}\n\\footnotesize"));
% f = strrep(f,"legend style={","legend style={font=\tiny,");
f = strrep(f,"mark size=0.5000pt","mark size=1.0000pt");

fid  = fopen(tikzpath,'w');
fprintf(fid,'%s',f);
fclose(fid);

end