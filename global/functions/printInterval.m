function printInterval(I)
% printInterval - prints an interval such that if one executes this command
% in the workspace, this intreval would be created
%
% Syntax:  
%    printInterval(I)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      12-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

fprintf('%s\n','interval( ...');

%print infimum
inf = infimum(I);
%first element
fprintf('[');
fprintf('%16.16f',inf(1));
fprintf('; ');
%other elements
if length(I)>2
    for i=2:(length(I)-1)
        fprintf('%16.16f',inf(i));
        fprintf('; ');
    end
else
    i = 1;
end
%last element
fprintf('%16.16f',inf(i+1));
fprintf('%s\n', '], ...');

%print supremum
sup = supremum(I);
%first element
fprintf('[');
fprintf('%16.16f',sup(1));
fprintf('; ');
if length(I)>2
    for i=2:(length(I)-1)
        fprintf('%16.16f',sup(i));
        fprintf('; ');
    end
else
    i = 1;
end
%close expression
fprintf('%16.16f',sup(i+1));
fprintf(']);');

%------------- END OF CODE --------------

