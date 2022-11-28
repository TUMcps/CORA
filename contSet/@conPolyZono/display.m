function display(cPZ)
% display - displays the properties of a conPolyZono object (center,
%    generator matrices, exponent matrix, constraint system) on the 
%    command window
%
% Syntax:  
%    display(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    ---
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

if isemptyobject(cPZ)
    
    dispEmptyObj(cPZ,inputname(1));

else

    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    % display center
    disp('c: ');
    disp(center(cPZ));
    
    % display generators
    G = cPZ.G;
    displayGenerators(G,DISPLAYDIM_MAX,'G');
    
    % display exponential matrix
    expMat = cPZ.expMat;
    displayGenerators(expMat,DISPLAYDIM_MAX,'expMat');
    
    % display constraint offset
    disp('b:');
    disp(cPZ.b);
    
    % display constraint generators
    A = cPZ.A;
    displayGenerators(A,DISPLAYDIM_MAX,'A');
    
    % display constraint exponential matrix
    expMat_ = cPZ.expMat_;
    displayGenerators(expMat_,DISPLAYDIM_MAX,'expMat_');
    
    % display independent generators
    Grest = cPZ.Grest;
    displayGenerators(Grest,DISPLAYDIM_MAX,'Grest');
    
    % display id
    disp('id:');
    disp(cPZ.id);

end

%------------- END OF CODE --------------