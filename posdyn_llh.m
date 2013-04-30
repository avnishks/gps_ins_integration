function [sys,x0,str,ts]=posdyn_llh(t,x,u,flag,lati,loni,alti)
%position dynamics - latitude, longitude & altitude.
%   
%   initial conditions (lati,loni,alti).
%   
   pdot = zeros(3,1);
  
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0         
    [sys,x0,str,ts] = mdlInitializeSizes(lati,loni,alti);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys = mdlDerivatives(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%%
  % Update and Terminate %
  %%%%%%%%%%%%%%%%%%%%%%%%
  case {2,9}
    sys = []; % do nothing

  %%%%%%%%%%
  % Output %
  %%%%%%%%%%
  case 3
    sys = mdlOutputs(t,x,u); 

  otherwise
    error(['unhandled flag = ',num2str(flag)]);
    
end
% end posdyn_llh
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes(lati,loni,alti)

sizes = simsizes;
sizes.NumContStates  = 3;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 3;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];
x0(1) = lati;
x0(2) = loni;
x0(3) = alti;
ts  = [0 0];   % sample time: [period, offset]

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Compute derivatives for continuous states.
%=============================================================================
%
function sys = mdlDerivatives(t,x,u)

   re = 20925604.;
   ecc = 8.1819191e-2;
   omes = 1.0 - ecc^2*sin(x(1))^2;
   R_n = re/sqrt(omes);
   R_m = R_n*(1.0-ecc^2)/(omes*sqrt(omes));

   pdot(1) = u(1)/(R_m+x(3));
   pdot(2) = u(2)/((R_n+x(3))*cos(x(1)));
   pdot(3) =-u(3);
   
   sys = pdot;

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)

   sys = x;

% end mdlOutputs