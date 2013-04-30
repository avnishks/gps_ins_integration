function [sys,x0,str,ts]=attdyn_DCM(t,x,u,flag,cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi)
%attitude dynamics - Direction Cosine Matrix (DCM).
%   
%   initial conditions (cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi).
%   

   cdot = zeros(9,1);
   
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0         
    [sys,x0,str,ts] = mdlInitializeSizes(cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi);

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
% end attdyn_DCM
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes(cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi)

sizes = simsizes;
sizes.NumContStates  = 9;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 9;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];
x0(1) = cxxi;
x0(2) = cxyi;
x0(4) = cyxi;
x0(5) = cyyi;
x0(7) = czxi;
x0(8) = czyi;
x0(3) = cyxi*czyi-cyyi*czxi;
x0(6) = cxyi*czxi-cxxi*czyi;
x0(9) = cxxi*cyyi-cxyi*cyxi;
ts  = [0 0];   % sample time: [period, offset]

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Compute derivatives for continuous states.
%=============================================================================
%
function sys = mdlDerivatives(t,x,u)

   cdot(1) = u(2)*x(7) - u(3)*x(4);
   cdot(2) = u(2)*x(8) - u(3)*x(5);
   cdot(3) = u(2)*x(9) - u(3)*x(6);
   cdot(4) = u(3)*x(1) - u(1)*x(7);
   cdot(5) = u(3)*x(2) - u(1)*x(8);
   cdot(6) = u(3)*x(3) - u(1)*x(9);
   cdot(7) = u(1)*x(4) - u(2)*x(1);
   cdot(8) = u(1)*x(5) - u(2)*x(2);
   cdot(9) = u(1)*x(6) - u(2)*x(3);
   
   sys = cdot;

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)

   x(3) = x(4)*x(8)-x(5)*x(7);
   x(6) = x(2)*x(7)-x(1)*x(8);
   x(9) = x(1)*x(5)-x(2)*x(4);

   sys = x;

% end mdlOutputs