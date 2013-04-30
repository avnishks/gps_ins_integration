function [sys,x0,str,ts] = gpsrcvr(t,x,u,flag,dt,p1i,p2i,p3i,p4i,p5i,p6i,p7i,p8i,p9i,p10i,p11i,p12i,vxni,vyni,vzni,rho_noise)
%GPS receiver Kalman filter.  
%
% GPS constelation's orbital parameters
   arg_of_perigee = zeros(18,1);
   right_ascension = zeros(18,1);
    
   x_pos_ECEF = zeros(4,1);
   y_pos_ECEF = zeros(4,1);
   z_pos_ECEF = zeros(4,1);
   rho_user = zeros(4,1);
   unit_los_x = zeros(4,1);
   unit_los_y = zeros(4,1);
   unit_los_z = zeros(4,1);
   
   State_Phi = zeros(12,12);
   Q_tilt_alt_clock = zeros(6);
   Q_pos = zeros(3,3);
   Q_vel = zeros(3,3);
   P_mat = zeros(12,12);
   pm = zeros(12,12);
   State_Vector = zeros(12,1);
   sv = zeros(12,1);
   resid_mat = zeros(4,4);
   K_mat = zeros(12,4);
   H_mat = zeros(4,12);
   R_mat = zeros(4,4);
   d_rho_sat = zeros(4,1);
   Y_mat = zeros(12,12);
   y = zeros(27,1);
  
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(dt,p1i,p2i,p3i,p4i,p5i,p6i,p7i,p8i,p9i,p10i,p11i,p12i,vxni,vyni,vzni);

  %%%%%%%%%%%%%%%
  % udate      %
  %%%%%%%%%%%%%%%
  case 2
    sys = mdlupdate(t,x,u,dt,rho_noise);
    
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys = mdlOutputs(t,x,u);
        
  %%%%%%%%%%%%%%%
  % udate      %
  %%%%%%%%%%%%%%%
  case 4
    sys = mdlGetTimeofNextVarHit(t,x,u,dt);
 
  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 1, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end
% end GPSrcvr
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(dt,p1i,p2i,p3i,p4i,p5i,p6i,p7i,p8i,p9i,p10i,p11i,p12i,vxni,vyni,vzni)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 12*(12+1)/2+12+3;
sizes.NumOutputs     = 12 + 12 + 3;
sizes.NumInputs      =24;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
str = [];

x0 = 0.0;
x0(1)  = p1i;
x0(13) = p2i;
x0(24) = p3i;
x0(34) = p4i;
x0(43) = p5i;
x0(51) = p6i;
x0(58) = p7i;
x0(64) = p8i;
x0(69) = p9i;
x0(73) = p10i;
x0(76) = p11i;
x0(78) = p12i;
x0(90) = 0.0;
x0(91) = vxni;
x0(92) = vyni;
x0(93) = vzni;
ts  = [dt 0];   % sample time: [period, offset]

% end mdlInitializeSizes
%
%=============================================================================
% mdludate
% Compute udates for discrete states.
%=============================================================================
%
function sys = mdlupdate(t,x,u,dt,rho_noise)

   Number_GPS_satellites = 18;
   irho_noise = rho_noise;
   
   arg_of_perigee(1)=    0.0;
   arg_of_perigee(2)=  120.0;  
   arg_of_perigee(3)=  240.0;
   arg_of_perigee(4)=   40.0;  
   arg_of_perigee(5)=  160.0;  
   arg_of_perigee(6)=  280.0;
   arg_of_perigee(7)=   80.0;  
   arg_of_perigee(8)=  200.0;  
   arg_of_perigee(9)=  320.0;
   arg_of_perigee(10)= 120.0;  
   arg_of_perigee(11)= 240.0;  
   arg_of_perigee(12)= 360.0;
   arg_of_perigee(13)= 160.0;  
   arg_of_perigee(14)= 280.0;  
   arg_of_perigee(15)=  40.0;
   arg_of_perigee(16)= 200.0;  
   arg_of_perigee(17)= 320.0;  
   arg_of_perigee(18)=  80.0;
   right_ascension(1)=    0.0;   
   right_ascension(2)=    0.0;   
   right_ascension(3)=    0.0; 
   right_ascension(4)=   60.0;   
   right_ascension(5)=   60.0;   
   right_ascension(6)=   60.0; 
   right_ascension(7)=  120.0;   
   right_ascension(8)=  120.0;   
   right_ascension(9)=  120.0; 
   right_ascension(10)= 180.0;   
   right_ascension(11)= 180.0;   
   right_ascension(12)= 180.0; 
   right_ascension(13)= 240.0;   
   right_ascension(14)= 240.0;   
   right_ascension(15)= 240.0; 
   right_ascension(16)= 300.0;   
   right_ascension(17)= 300.0;   
   right_ascension(18)= 300.0;

% constants
   pi = 3.141592654;
   deg_to_rad = 0.01745329252;
   rad_to_deg = 57.29577951;
   meters_to_feet = 3.2808398;
   R_equator = 20925604.0;
   ellipticity = 3.35281066474E-03;
   eccentricity = 8.18191908426E-02;
   omega_earth = 7.292115E-05;

% GPS constants
   grav_mu_WGS84 = 1.4076443E+16;
   semi_major = 87142134.0;
   GPS_eccentricity = 0.005;
   orbit_incline = 55.0;
   
   meters_to_feet = 3.2808398;
   NumberStates = 12;
   tau_altimeter = 1000.0;
   
   elp = 3.3528107e-3;
   g0 = 32.087641;
   g1 = -2.01364;
   g2 = 5.31658e-3;
   
      sim_time = t;
      user_lat = u(1);
      user_lon = u(2);
      user_alt = u(3);
      user_v_n = u(4);
      user_v_e = u(5);
      user_v_u = u(6);
      user_alpha = u(16);

      k = 0;
      for i=1:NumberStates
          for j=i:NumberStates
             k = k+1;
             P_mat(i,j) = x(k);
             P_mat(j,i) = P_mat(i,j);
         end
         State_Vector(i,1) = x(i+78);
      end

      s_user_lat = sin ( user_lat*deg_to_rad );
      c_user_lat = cos ( user_lat*deg_to_rad ); 
      s_user_lon = sin ( user_lon*deg_to_rad ); 
      c_user_lon = cos ( user_lon*deg_to_rad );
      s_user_alpha = sin ( user_alpha );
      c_user_alpha = cos ( user_alpha );

% form the local-level to ECEF direction cosine matrix

      CLExx = - c_user_alpha*s_user_lon - s_user_alpha*c_user_lon*s_user_lat;
      CLExy =   s_user_alpha*s_user_lon - c_user_alpha*c_user_lon*s_user_lat;
      CLExz =   c_user_lon*c_user_lat;
      CLEyx =   c_user_alpha*c_user_lon - s_user_alpha*s_user_lon*s_user_lat;
      CLEyy = - s_user_alpha*c_user_lon - c_user_alpha*s_user_lon*s_user_lat;
      CLEyz =   s_user_lon*c_user_lat;
      CLEzx =   s_user_alpha*c_user_lat;
      CLEzy =   c_user_alpha*c_user_lat;
      CLEzz =   s_user_lat;

% user position in ECEF

      sqr_one_f = ( 1.0 - ellipticity )^2;
      denom = sqrt ( c_user_lat^2 + sqr_one_f*s_user_lat^2 );
      R = R_equator/denom + user_alt;
      S = R_equator*sqr_one_f/denom + user_alt;

      user_x_ECEF = R*c_user_lat*c_user_lon; 
      user_y_ECEF = R*c_user_lat*s_user_lon; 
      user_z_ECEF = S*s_user_lat;

      user_x_vel_nav =  c_user_alpha*user_v_n - s_user_alpha*user_v_e;
      user_y_vel_nav = -s_user_alpha*user_v_n - c_user_alpha*user_v_e;
      user_z_vel_nav =                                                  user_v_u;

      vxd = ( user_x_vel_nav - x(91) )/dt;
      vyd = ( user_y_vel_nav - x(92) )/dt;
      vzd = ( user_z_vel_nav - x(93) )/dt;
      x(91) = user_x_vel_nav;
      x(92) = user_y_vel_nav;
      x(93) = user_z_vel_nav;

      ovrRy = ( 1.0 - user_alt/R_equator + elp*(2.0*CLExy^2-CLExz^2) )/R_equator;
      ovrRx = ( 1.0 - user_alt/R_equator + elp*(2.0*CLExx^2-CLExz^2) )/R_equator;
      ovrT  =   2.0*elp*CLExx*CLExy/R_equator;
      rhox  = -user_y_vel_nav*ovrRy - user_x_vel_nav*ovrT;
      rhoy  =  user_x_vel_nav*ovrRx + user_y_vel_nav*ovrT;
%   
      omx = CLExx*omega_earth;
      omy = CLExy*omega_earth;
      omz = CLExz*omega_earth;
%   
      cnx = (rhoy+2.0*omy)*user_z_vel_nav - (     2.0*omz)*user_y_vel_nav;
      cny = (     2.0*omz)*user_x_vel_nav - (rhox+2.0*omx)*user_z_vel_nav;
      cnz = (rhox+2.0*omx)*user_y_vel_nav - (rhoy+2.0*omy)*user_x_vel_nav;
%
      gnz = -g0*(1.0+g1*(user_alt/R_equator)+g2*CLExz^2);
%   
      f_x_nav = vxd + cnx;
      f_y_nav = vyd + cny;
      f_z_nav = vzd + cnz - gnz;

%  form state transition matrix

      State_Phi(1,1) = 1.0;
      State_Phi(1,3) = -rhoy*dt;
      State_Phi(1,4) = dt;
      State_Phi(2,2) = 1.0;
      State_Phi(2,3) = rhox*dt;
      State_Phi(2,5) = dt;
      State_Phi(3,1) = -State_Phi(1,3);
      State_Phi(3,2) = -State_Phi(2,3);
      State_Phi(3,3) = 1.0;
      State_Phi(3,6) = dt;
      State_Phi(4,1) = gnz*dt/R_equator;
      State_Phi(4,4) = 1.0;
      State_Phi(4,5) = 2.0*omz*dt;
      State_Phi(4,6) = -(rhoy+2.0*omy)*dt;
      State_Phi(4,8) = -f_z_nav*dt;
      State_Phi(4,9) = f_y_nav*dt;
      State_Phi(5,2) = gnz*dt/R_equator;
      State_Phi(5,4) = -State_Phi(4,5);
      State_Phi(5,5) = 1.0;
      State_Phi(5,6) = (rhox+2.0*omx)*dt;
      State_Phi(5,7) = -State_Phi(4,8);
      State_Phi(5,9) = -f_x_nav*dt;
      State_Phi(6,3) = -2.0*gnz*dt/R_equator;
      State_Phi(6,4) = -State_Phi(4,6);
      State_Phi(6,5) = -State_Phi(5,6);
      State_Phi(6,6) = 1.0;
      State_Phi(6,7) = -State_Phi(4,9);
      State_Phi(6,8) = -State_Phi(5,9);
      State_Phi(7,7) = 1.0;
      State_Phi(7,8) = omz*dt;
      State_Phi(7,9) = -(rhoy+omy)*dt;
      State_Phi(8,7) = -State_Phi(7,8);
      State_Phi(8,8) = 1.0;
      State_Phi(8,9) = (rhox+omx)*dt;
      State_Phi(9,7) = -State_Phi(7,9);
      State_Phi(9,8) = -State_Phi(8,9);
      State_Phi(9,9) = 1.0;
      State_Phi(10,10) = exp( -dt/tau_altimeter );
      State_Phi(11,11) = 1.0;
      State_Phi(11,12) = dt;
      State_Phi(12,12) = 1.0;

%  state vector propagation

      sv = State_Phi*State_Vector;
      State_Vector = sv;
      
%  compute process noise matrix components

%  process noise

      Q_LL_ph = 1.0*( meters_to_feet^2 );
      Q_LL_pv = 1.0*( meters_to_feet^2 );
      Q_LL_vh = 1.0E-3*( meters_to_feet^2 );
      Q_LL_vv = 1.0E-2*( meters_to_feet^2 );
      Q_tilt_alt_clock(1) = 1.0E-12;
      Q_tilt_alt_clock(2) = 1.0E-12;
      Q_tilt_alt_clock(3) = 1.0E-12;
      Q_tilt_alt_clock(4) = 10.0*( meters_to_feet^2 );
      Q_tilt_alt_clock(5) = 0.25*( meters_to_feet^2 );
      Q_tilt_alt_clock(6) = 2.0E-3*( meters_to_feet^2);
       
      Q_pos(1,1) = Q_LL_ph*dt;
      Q_pos(2,2) = Q_LL_ph*dt;
      Q_pos(3,3) = Q_LL_pv*dt;
      Q_vel(1,1) = Q_LL_vh*dt;
      Q_vel(2,2) = Q_LL_vh*dt;
      Q_vel(3,3) = Q_LL_vv*dt;

%  propagate error covariance matrix

      pm = State_Phi*P_mat*State_Phi';
      P_mat = pm;

      for jstate = 1:3
         for istate = 1:3
            P_mat(istate,jstate) = P_mat(istate,jstate) + Q_pos(istate,jstate);
            P_mat(istate+3,jstate+3) = P_mat(istate+3,jstate+3) + Q_vel(istate,jstate);
         end
      end

      for istate = 1:6
         P_mat(istate+6,istate+6) = P_mat(istate+6,istate+6) + Q_tilt_alt_clock(istate)*dt;
      end

% begin satellite orbit position computations

      mean_motion = sqrt ( grav_mu_WGS84/ (semi_major*semi_major*semi_major) );
      mean_anomaly = mean_motion*sim_time;

% initial guess of eccentric anomaly for ineration

      eccen_anomaly = mean_anomaly + GPS_eccentricity*sin( mean_anomaly );

% iteration on eccentric anomaly

      for i_iter=1:2
         s_eccen = sin ( eccen_anomaly );
         c_eccen = cos ( eccen_anomaly );
         eccen_anomaly = ( GPS_eccentricity*(s_eccen-eccen_anomaly*c_eccen)+mean_anomaly )/(1.0-GPS_eccentricity*c_eccen);
      end

      s_eccen = sin ( eccen_anomaly );
      c_eccen = cos ( eccen_anomaly );
      c_true_anomaly = (c_eccen-GPS_eccentricity)/(1.0-GPS_eccentricity*c_eccen);
      one_eccens = 1.0 - GPS_eccentricity*GPS_eccentricity;
      s_true_anomaly = sqrt( one_eccens )*s_eccen/(1.0-GPS_eccentricity*c_eccen);

      true_anomaly = atan( s_true_anomaly/c_true_anomaly );
      if ( s_true_anomaly>0.0 ) & ( c_true_anomaly<0.0 )
         true_anomaly = pi - atan( s_true_anomaly/(-c_true_anomaly) );
      end
      if ( s_true_anomaly<0.0 ) & ( c_true_anomaly<0.0 ) 
         true_anomaly = atan( -s_true_anomaly/(-c_true_anomaly) ) - pi;
      end      
      if ( s_true_anomaly<0.0 ) & ( c_true_anomaly>0.0 )
         true_anomaly = -atan( -s_true_anomaly/c_true_anomaly );
      end
      
      orbit_radius = semi_major*( 1.0 - GPS_eccentricity*c_eccen );

% compute relative range based on 4 selected satellites

      for j_sat=1:4
          
         i_sat = u(j_sat+16);
          
         arg_latitude = true_anomaly + arg_of_perigee(i_sat)*deg_to_rad;

         s_arg_lat = sin( arg_latitude );
         c_arg_lat = cos( arg_latitude ); 

         x_pos_orbit = orbit_radius*c_arg_lat;
         y_pos_orbit = orbit_radius*s_arg_lat; 

         lon_ascend_node = right_ascension(i_sat)*deg_to_rad - omega_earth*sim_time;

         s_lon_ascend = sin( lon_ascend_node ); 
         c_lon_ascend = cos( lon_ascend_node );

         s_incline = sin( orbit_incline*deg_to_rad );
         c_incline = cos( orbit_incline*deg_to_rad );

         x_pos_ECEF(j_sat) = x_pos_orbit*c_lon_ascend - y_pos_orbit*c_incline*s_lon_ascend;
         y_pos_ECEF(j_sat) = x_pos_orbit*s_lon_ascend + y_pos_orbit*c_incline*c_lon_ascend; 
         z_pos_ECEF(j_sat) = y_pos_orbit*s_incline;

         los_x_ECEF = (x_pos_ECEF(j_sat) ) - user_x_ECEF; 
         los_y_ECEF = (y_pos_ECEF(j_sat) ) - user_y_ECEF; 
         los_z_ECEF = (z_pos_ECEF(j_sat) ) - user_z_ECEF; 

         rho_user(j_sat) = sqrt ( los_x_ECEF^2 + los_y_ECEF^2 + los_z_ECEF^2 );

         los_x_nav = CLExx*los_x_ECEF + CLEyx*los_y_ECEF + CLEzx*los_z_ECEF;
         los_y_nav = CLExy*los_x_ECEF + CLEyy*los_y_ECEF + CLEzy*los_z_ECEF;
         los_z_nav = CLExz*los_x_ECEF + CLEyz*los_y_ECEF + CLEzz*los_z_ECEF;
         
         measr_rho          = u(j_sat+20)  + irho_noise(floor(t) + 1,j_sat) ;
         d_rho_sat(j_sat,1) = measr_rho - (rho_user(j_sat) );

         unit_los_x(j_sat) = los_x_nav/rho_user(j_sat);
         unit_los_y(j_sat) = los_y_nav/rho_user(j_sat);
         unit_los_z(j_sat) = los_z_nav/rho_user(j_sat);
         
         H_mat(j_sat,1) = unit_los_x(j_sat);
         H_mat(j_sat,2) = unit_los_y(j_sat);
         H_mat(j_sat,3) = unit_los_z(j_sat);
         H_mat(j_sat,11) = 1.0;
         H_mat(j_sat,12) = 0.0;
        
      end
 
%  measurement updates - psuedo range

      R_var_rho = 10.0*( meters_to_feet^2 );
      R_mat(1,1) = R_var_rho;
      R_mat(2,2) = R_var_rho;
      R_mat(3,3) = R_var_rho;
      R_mat(4,4) = R_var_rho;
      Y_mat = eye(NumberStates);

      resid_mat = H_mat*P_mat*H_mat' + R_mat;
      K_mat = P_mat*H_mat'*inv(resid_mat);
      pm = (Y_mat - K_mat*H_mat)*P_mat;
      P_mat = pm;
      sv = State_Vector + K_mat*(d_rho_sat - H_mat*State_Vector);
      State_Vector = sv;

      k = 0;
      for i=1:NumberStates
         for j=i:NumberStates
             k = k+1;
             x(k) = P_mat(i,j);
         end
         x(i+78) = State_Vector(i,1);
      end
   
   sys = x;

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys = mdlOutputs(t,x,u)

   R_equator = 20925604.0;

      user_v_n = u(4);
      user_v_e = u(5);
      user_v_u = u(6);
      user_alpha = u(16);
      
      s_user_alpha = sin ( user_alpha );
      c_user_alpha = cos ( user_alpha );
      
      user_lat = u(1);
      user_lon = u(2);
      user_alt = u(3);
      deg_to_rad = (pi/180);
      s_user_lat = sin ( user_lat*deg_to_rad );
      c_user_lat = cos ( user_lat*deg_to_rad ); 
      s_user_lon = sin ( user_lon*deg_to_rad ); 
      c_user_lon = cos ( user_lon*deg_to_rad );
      CLExx = - c_user_alpha*s_user_lon - s_user_alpha*c_user_lon*s_user_lat;
      CLExy =   s_user_alpha*s_user_lon - c_user_alpha*c_user_lon*s_user_lat;
      CLExz =   c_user_lon*c_user_lat;
      CLEyx =   c_user_alpha*c_user_lon - s_user_alpha*s_user_lon*s_user_lat;
      CLEyy = - s_user_alpha*c_user_lon - c_user_alpha*s_user_lon*s_user_lat;
      CLEyz =   s_user_lon*c_user_lat;
      CLEzx =   s_user_alpha*c_user_lat;
      CLEzy =   c_user_alpha*c_user_lat;
      CLEzz =   s_user_lat;
      

      user_x_vel_nav =  c_user_alpha*user_v_n - s_user_alpha*user_v_e;
      user_y_vel_nav = -s_user_alpha*user_v_n - c_user_alpha*user_v_e;
      user_z_vel_nav =                                                  user_v_u;
      
      dvxn = x(82) - user_z_vel_nav*x(79)/R_equator;
      dvyn = x(83)                                  - user_z_vel_nav*x(80)/R_equator;
      dvzn = x(84) + user_y_vel_nav*x(80)/R_equator + user_x_vel_nav*x(79)/R_equator;
      
      y(1) =  c_user_alpha*dvxn - s_user_alpha*dvyn;
      y(2) =  s_user_alpha*dvxn + c_user_alpha*dvyn;
      y(3) =                                         -dvzn;

      rcvr_psi_L_x = x(85);
      rcvr_psi_L_y = x(86);
      rcvr_psi_L_z = x(87);

      y(4) = rcvr_psi_L_x - x(80)/R_equator;
      y(5) = rcvr_psi_L_y + x(79)/R_equator;
      y(6) =-rcvr_psi_L_z;
      
      y(7) = CLExx*x(79) + CLExy*x(80) + CLExz*x(81);
      y(8) = CLEyx*x(79) + CLEyy*x(80) + CLEyz*x(81);
      y(9) = CLEzx*x(79) + CLEzy*x(80) + CLEzz*x(81);
      y(10) = x(88);
      y(11) = x(89);
      y(12) = x(90);
      
      y(13) = x(1);
      y(14) = x(13);
      y(15) = x(24);
      y(16) = x(34);
      y(17) = x(43);
      y(18) = x(51);
      y(19) = x(58);
      y(20) = x(64);
      y(21) = x(69);
      y(22) = x(73);
      y(23) = x(76);
      y(24) = x(78);
      
      y(25) = x(79);
      y(26) = x(80);
      y(27) = x(81);

  sys = y;
  
% end mdlOutputs
%
%=============================================================================
% mdlGetTimeofNextVarHit
%=============================================================================
function sys = mdlGetTimeofNextVarHit(t,x,u,dt)

   sys = t + dt;
  
% end mdlGetTimeofNextVarHit
