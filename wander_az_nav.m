function [sys,x0,str,ts]=wander_az_nav(t,x,u,flag,ai,bi,ci,di,vxi,vyi,vzi,cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi,hi)
%wander azimuth INU dynamics.
%   
%   initial conditions (ai,bi,ci,di,vxi,vyi,vzi,cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi,hi).
%   
   xdot = zeros(14,1);
   y = zeros(16,1);
   
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0         
    [sys,x0,str,ts] = mdlInitializeSizes(ai,bi,ci,di,vxi,vyi,vzi,cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi,hi);

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
function [sys,x0,str,ts] = mdlInitializeSizes(ai,bi,ci,di,vxi,vyi,vzi,cxxi,cxyi,cxzi,cyxi,cyyi,cyzi,czxi,czyi,czzi,hi)

sizes = simsizes;
sizes.NumContStates  = 14;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 16;
sizes.NumInputs      = 7;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];

x0(1) = ai;
x0(2) = bi;
x0(3) = ci;
x0(4) = di;

x0(5) = vxi;
x0(6) = vyi;
x0(7) = vzi;

x0(8)  = cxxi;
x0(9)  = cxyi;
x0(10) = cxzi;
x0(11) = cyxi;
x0(12) = cyyi;
x0(13) = cyzi;

x0(14) = hi;

ts  = [0 0];   % sample time: [period, offset]

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Compute derivatives for continuous states.
%=============================================================================
%
function sys = mdlDerivatives(t,x,u,anx,any,anz)
%
% gyro inputs
   wbx = u(1);
   wby = u(2);
   wbz = u(3);
% accelerometer inputs   
   abx = u(4);
   aby = u(5);
   abz = u(6);
% altitude input
   alt = u(7);
   
% quaternion variables   
   a   = x(1);
   b   = x(2);
   c   = x(3);
   d   = x(4);
% velocity variables
   vx  = x(5);
   vy  = x(6);
   vz  = x(7);
% CNE DCM variables
   cnxx = x(8);
   cnxy = x(9);
   cnxz = x(10);
   cnyx = x(11);
   cnyy = x(12);
   cnyz = x(13);
% altitude
   h   = x(14);
%   
% constants
   re  = 20925604.;
   elp = 3.3528107e-3;
   ome = 7.292115e-5;
   g0 = 32.087641;
   g1 = -2.01364;
   g2 = 5.31658e-3;
   k1 = 1.0e-4;
   k2 = 1.0e-6;
%
   ovrRy = ( 1.0 - h/re + elp*(2.0*cnxy^2-cnxz^2) )/re;
   ovrRx = ( 1.0 - h/re + elp*(2.0*cnxx^2-cnxz^2) )/re;
   ovrT  =   2.0*elp*cnxx*cnxy/re;
   rhox  = -vy*ovrRy - vx*ovrT;
   rhoy  =  vx*ovrRx + vy*ovrT;
%   
   omx = cnxx*ome;
   omy = cnxy*ome;
   omz = cnxz*ome;
%   
% normalize quaternions
   rss = sqrt(a^2 + b^2 + c^2 + d^2);
   a = a/rss;
   b = b/rss;
   c = c/rss;
   d = d/rss;
%
   cbxx =  (a*a + d*d - b*b - c*c);
   cbyx = -2.0*(a*b + c*d);
   cbzx = -2.0*(a*c - b*d);
   cbxy =  2.0*(a*b - c*d);
   cbyy = -(b*b + d*d - a*a - c*c);
   cbzy = -2.0*(b*c + a*d);
   cbxz = cbyx*cbzy-cbyy*cbzx;
   cbyz = cbxy*cbzx-cbxx*cbzy;
   cbzz = cbxx*cbyy-cbxy*cbyx;
%
   wcx = wbx - (cbxx*omx+cbyx*omy+cbzx*omz) - (cbxx*rhox+cbyx*rhoy);
   wcy = wby - (cbxy*omx+cbyy*omy+cbzy*omz) - (cbxy*rhox+cbyy*rhoy);
   wcz = wbz - (cbxz*omx+cbyz*omy+cbzz*omz) - (cbxz*rhox+cbyz*rhoy);
%
   adot = 0.5*(         wcz*b - wcy*c + wcx*d );
   bdot = 0.5*(-wcz*a         + wcx*c + wcy*d );
   cdot = 0.5*( wcy*a - wcx*b         + wcz*d );
   ddot = 0.5*(-wcx*a - wcy*b - wcz*c         );
%
   cnx = (rhoy+2.0*omy)*vz - (     2.0*omz)*vy;
   cny = (     2.0*omz)*vx - (rhox+2.0*omx)*vz;
   cnz = (rhox+2.0*omx)*vy - (rhoy+2.0*omy)*vx;
%
   anx = cbxx*abx + cbxy*aby + cbxz*abz;
   any = cbyx*abx + cbyy*aby + cbyz*abz;
   anz = cbzx*abx + cbzy*aby + cbzz*abz;
%   
   gnz = -g0*(1.0+g1*(h/re)+g2*cnxz^2);
%   
   vxdot = anx - cnx;
   vydot = any - cny;
   vzdot = anz - cnz + gnz + k2*(alt-h);
%
   hdot = vz + k1*(alt-h);
%
   cxxdot =           - cnxz*rhoy;
   cxydot = cnxz*rhox;
   cxzdot = cnxx*rhoy - cnxy*rhox;
   cyxdot =           - cnyz*rhoy;
   cyydot = cnyz*rhox;
   cyzdot = cnyx*rhoy - cnyy*rhox;
%   
   xdot(1) = adot;
   xdot(2) = bdot;
   xdot(3) = cdot;
   xdot(4) = ddot;
   xdot(5) = vxdot;
   xdot(6) = vydot;
   xdot(7) = vzdot;
   xdot(8)  = cxxdot;
   xdot(9)  = cxydot;
   xdot(10) = cxzdot;
   xdot(11) = cyxdot;
   xdot(12) = cyydot;
   xdot(13) = cyzdot;
   xdot(14) = hdot;
%   
   sys = xdot;

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
%
function sys = mdlOutputs(t,x,u)

   rad_to_deg = 57.29577951;

   y(1) =  rad_to_deg*asin(x(10));
   cnzz = x(8)*x(12) - x(11)*x(9);
   y(2) =  rad_to_deg*atan(-x(13)/cnzz);
   y(3) =  x(14);
   alpha = atan(-x(9)/x(8));
   y(4) =  x(5)*cos(alpha)-x(6)*sin(alpha);
   y(5) = -x(5)*sin(alpha)-x(6)*cos(alpha);
   y(6) =  x(7);
   a   = x(1);
   b   = x(2);
   c   = x(3);
   d   = x(4);
   cbxx =  (a*a + d*d - b*b - c*c);
   cbyx = -2.0*(a*b + c*d);
   cbzx = -2.0*(a*c - b*d);
   cbxy =  2.0*(a*b - c*d);
   cbyy = -(b*b + d*d - a*a - c*c);
   cbzy = -2.0*(b*c + a*d);
   cbxz = cbyx*cbzy-cbyy*cbzx;
   cbyz = cbxy*cbzx-cbxx*cbzy;
   cbzz = cbxx*cbyy-cbxy*cbyx;
   y(7) = cbxx;
   y(8) = cbxy;
   y(9) = cbxz;
   y(10) = cbyx;
   y(11) = cbyy;
   y(12) = cbyz;
   y(13) = cbzx;
   y(14) = cbzy;
   y(15) = cbzz;
   y(16)= alpha;

   sys = y;

% end mdlOutputs