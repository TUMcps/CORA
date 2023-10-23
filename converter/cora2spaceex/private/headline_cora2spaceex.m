function [docNode,component] = headline_cora2spaceex()
% headline_cora2spaceex - convert CORA models to SpaceEx models
%
% Syntax:
%    [docNode,component] = headline_cora2spaceex()
%
% Inputs:
%    -
%
% Outputs:
%    docNode   - xml-file that contains the SpaceEx-model
%    component - 
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Create the document node and root element (spaceex)
docNode = com.mathworks.xml.XMLUtils.createDocument('sspaceex');

%Identify the root element (spaceex) and set the version and math attribute
spaceex = docNode.getDocumentElement;
spaceex.setAttribute('version','2.0');
spaceex.setAttribute('math','spaceex');

%Add the element node (component), for the product page and set the id
%attribute. This element node is the child of the root element(spaceex)
component = docNode.createElement('component');
component.setAttribute('id','model');
spaceex.appendChild(component);

end

% ------------------------------ END OF CODE ------------------------------
