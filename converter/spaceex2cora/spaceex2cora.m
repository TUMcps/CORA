function spaceex2cora(xmlFile,varargin)
% spaceex2cora - convert SpaceEx models to CORA models 
%
% Syntax:  
%    spaceex2cora(xmlFile)
%    spaceex2cora(xmlFile,save2pha,rootID,outputName,outputDir)
%
% Inputs:
%    xmlFile - path to the xml-file that contains the SpaceEx-model
%    save2pha - conversion to parallel hybrid automaton (default: 1)
%               conversion to flat hybrid automaton (0)
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

    % decision whether parallel HA or flat HA
    save2pha = 1; % default: parallel HA
    if nargin >= 2
        save2pha = varargin{1};
    end
    % parse the SpaceEX-model
    if nargin >= 4
        data = SX2structHA(xmlFile,save2pha,varargin{2},varargin{3});
    elseif nargin >= 3
        data = SX2structHA(xmlFile,save2pha,varargin{2});
    else
        data = SX2structHA(xmlFile,save2pha);
    end
    
    % generate the CORA model files
    if nargin >= 4
       if nargin >= 5
          StructHA2file(data,varargin{3},varargin{4});
       else
          StructHA2file(data,varargin{3});
       end
    else
       % default name == xml-file name
       [~,name,~] = fileparts(xmlFile); 
       
       % generate files
       StructHA2file(data,name);
    end
    
%------------- END OF CODE --------------