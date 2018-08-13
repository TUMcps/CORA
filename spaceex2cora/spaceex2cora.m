function spaceex2cora(xmlFile,varargin)
% spaceex2cora - convert SpaceEx models to CORA models 
%
% Syntax:  
%    spaceex2cora(xmlFile)
%    spaceex2cora(xmlFile,rootID,outputName,outputDir)
%
% Inputs:
%    xmlFile - path to the xml-file that contains the SpaceEx-model
%    rootID -  ID of SpaceEx component to be used as root component
%              (specified as a string)
%    outputName - name of the generated CORA model (specified as string)
%    outputDir - path to the desired output directory where all generated
%                files are stored
%
% Outputs:
%   ---
%
% Example: 
%    spaceex2cora('build_48.xml','sys');
%    cd([coraroot,'/models/SpaceExConverted/build_48']);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      03-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse the SpaceEX-model
    if nargin >= 2
       data = SX2structHA(xmlFile,varargin{1});
    else
       data = SX2structHA(xmlFile);
    end
    
    % genearate the CORA model files
    if nargin >= 3
       if nargin >= 4
          StructHA2file(data,varargin{2},varargin{3});
       else
          StructHA2file(data,varargin{2});
       end
    else
       % default name == xml-file name
       [~,name,~] = fileparts(xmlFile); 
       
       % generate files
       StructHA2file(data,name);
    end
    
%------------- END OF CODE --------------