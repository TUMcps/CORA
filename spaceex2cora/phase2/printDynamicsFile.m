function printDynamicsFile(path,name,dynamics)

auxpath = [path '/auxiliary'];
if ~exist(auxpath,'dir')
    mkdir(auxpath);
end

addpath(auxpath);

fname = strcat(auxpath,'/',name,'.m');
file = fopen(fname,'w');

if file<0
    error('could not open file "%s"',fname);
end

% write file contents
Str = "function [dx]=" + name + "(t,x,u)" + newline + newline + dynamics;
fwrite(file,Str);
fclose(file);

% ensure matlab detects new function
rehash path;

end