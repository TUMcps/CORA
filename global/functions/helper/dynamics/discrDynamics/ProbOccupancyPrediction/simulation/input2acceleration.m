function acceleration = input2acceleration(input,velocity,type)
% input2acceleration - transforms the input (-1,1) to a physical acceleration. 
%
% Syntax:
%    acceleration = input2acceleration(input,velocity,type)
%
% Inputs:
%    input - input of the vehicle ranging from (-1,1)
%    velocity - velocity of the vehicle
%
% Outputs:
%    acceleration - acceleration value of the car
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-June-2008
% Last update:   15-June-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(type,'car')
    vlong=7.32;
    af=9.1;
    k=66.6;
elseif strcmp(type,'bicycle')
    vlong=0.2;
    af=4;
    k=0.75;
end   


%preallocation
acceleration=zeros(size(input));

%input>=0 =>  a=af*u if v<=vlong   a=u*k/v if v>vlong
acceleration(velocity<=vlong & input>=0)=af.*input(velocity<=vlong &  input>=0);
acceleration(velocity>vlong & input>=0)=input(velocity>vlong & input>=0).*k./velocity(velocity>vlong & input>=0);

%input<0 =>     a=af*u if v>0      a=0 if v<=0 
acceleration(velocity>0 & input<0)=af.*input(velocity>0 & input<0);

% if vectorization is not necessary you can comment the lines above and use
% the following instead
%
% if (input>=0)
%     if(v<=vlong)
%         acceleration=af*input;
%     else
%         acceleration=k/v*input;
%     end
% else
%     if(v>0)
%         acceleration=af*input;
%     else
%         acceleration=0;
%     end 
% end
    
% ------------------------------ END OF CODE ------------------------------
