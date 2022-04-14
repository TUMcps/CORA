function display(cPZ)
% display - displays a conPolyZono object to the console
%
% Syntax:  
%    display(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    -
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    expMat = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    expMat_ = [0 1 2; 1 0 0; 0 1 0];
%
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/display

% Author:        Niklas Kochdumper
% Written:       19-January-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    % display center
    disp('c: ');
    disp(center(cPZ));
    
    % display generators
    G = cPZ.G;
    maxGens = 10;
    displayGenerators(G,maxGens,'G');
    
    % display exponential matrix
    expMat = cPZ.expMat;
    displayGenerators(expMat,maxGens,'expMat');
    
    % display constraint offset
    disp('b:');
    disp(cPZ.b);
    
    % display constraint generators
    A = cPZ.A;
    displayGenerators(A,maxGens,'A');
    
    % display constraint exponential matrix
    expMat_ = cPZ.expMat_;
    displayGenerators(expMat_,maxGens,'expMat_');
    
    % display independent generators
    Grest = cPZ.Grest;
    displayGenerators(Grest,maxGens,'Grest');
    
    % display id
    disp('id:');
    disp(cPZ.id);
    
end

%------------- END OF CODE --------------