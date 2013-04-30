function [sys,x0,str,ts] = gpssat(t,x,u,flag)
%GPS satellites' visibility and relative range.  
%
% GPS constelation's orbital parameters

   arg_of_perigee = zeros(18,1);
   right_ascension = zeros(18,1);
    
   x_pos_ECEF = zeros(18,1);
   y_pos_ECEF = zeros(18,1);
   z_pos_ECEF = zeros(18,1);
   rho_user = zeros(18,1);
   unit_los_x = zeros(18,1);
   unit_los_y = zeros(18,1);
   unit_los_z = zeros(18,1);
   elev_angle = zeros(18,1);
   unit_LL_x = zeros(18,1);
   unit_LL_y = zeros(18,1);
   unit_LL_z = zeros(18,1);
   
   i_sat_vis = zeros(18,1);
   G_matrix = zeros(4,4);
   GDOP_matrix = zeros(4,4);
      
   y = zeros(14);

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 1, 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end
% end GPSrange
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     =14;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [-1 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

   Number_GPS_satellites = 18;
   
   arg_of_perigee(1)=    0.0;
   arg_of_perigee(2)=  120.0;  
   arg_of_perigee(3)=  240.0;
   arg_of_perigee(4)=   40.0;  
   arg_of_perigee(5)=  160.0;  
   arg_of_perigee(6)=  280.0;
   arg_of_perigee(7)=   80.0;  
   arg_of_perigee(8)=  200.0;  
   arg_of_perigee(9)=  320.0;
   arg_of_perigee(10)=  120.0;  
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
   right_ascension(10)=  180.0;   
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

      sim_time = t;
      user_lat = u(1)*rad_to_deg;
      user_long = u(2)*rad_to_deg;
      user_alt = u(3);

      s_user_lat = sin ( user_lat*deg_to_rad );
      c_user_lat = cos ( user_lat*deg_to_rad ); 
      s_user_long = sin ( user_long*deg_to_rad ); 
      c_user_long = cos ( user_long*deg_to_rad );

% form the local-level to ECEF direction cosine matrix

      CLExx = - s_user_lat*c_user_long;
      CLExy = - s_user_long;
      CLExz = - c_user_lat*c_user_long;
      CLEyx = - s_user_lat*s_user_long; 
      CLEyy =   c_user_long;
      CLEyz = - c_user_lat*s_user_long;
      CLEzx =   c_user_lat;
      CLEzy =   0.0;
      CLEzz = - s_user_lat;

% user position in ECEF

      sqr_one_f = ( 1.0 - ellipticity )*( 1.0 - ellipticity );
      denom = sqrt ( (c_user_lat*c_user_lat) + sqr_one_f*(s_user_lat*s_user_lat) );
      R = R_equator/denom + user_alt;
      S = R_equator*sqr_one_f/denom + user_alt;

      user_x_ECEF = R*c_user_lat*c_user_long; 
      user_y_ECEF = R*c_user_lat*s_user_long; 
      user_z_ECEF = S*s_user_lat;

      R_user_ECEF = sqrt( user_x_ECEF^2 + user_y_ECEF^2 + user_z_ECEF^2 );

% local vertical at user's position in ECEF

      k2 = ( 1.0 - ellipticity )*( 1.0 - ellipticity );
      e2 = 1.0 - k2; 
      w_o = sqrt ( user_x_ECEF*user_x_ECEF + user_y_ECEF*user_y_ECEF );
      if  w_o >= 1.0
         R_N = R_equator;
         H = 0.0;
         Q_o = user_z_ECEF / w_o;
            for i_iter=1:2
               TAN1 = Q_o*( R_N + H )/( R_N*k2 + H );
               SIN2 = (TAN1*TAN1)/( 1.0 + TAN1*TAN1 ); 
               R_N = R_equator/sqrt( 1.0 - e2*SIN2 );
               H = w_o*sqrt( 1.0 + TAN1*TAN1 ) - R_N;
            end
      else
         R_N = R_equator*( 1.0 - ellipticity )/k2;
         H = abs( user_z_ECEF ) - R_N*k2;
      end

      up_user_x = user_x_ECEF/( R_N + H );
      up_user_y = user_y_ECEF/( R_N + H );
      up_user_z = user_z_ECEF/( R_N*k2 + H );

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

% identify candidate satellites based on elevation angle and distance */

      
      i_temp = 0;
      for i_sat=1:Number_GPS_satellites
         arg_latitude = true_anomaly + arg_of_perigee(i_sat)*deg_to_rad;

         s_arg_lat = sin( arg_latitude );
         c_arg_lat = cos( arg_latitude ); 

         x_pos_orbit = orbit_radius*c_arg_lat;
         y_pos_orbit = orbit_radius*s_arg_lat; 

         long_ascend_node = right_ascension(i_sat)*deg_to_rad - omega_earth*sim_time;

         s_long_ascend = sin( long_ascend_node ); 
         c_long_ascend = cos( long_ascend_node );

         s_incline = sin( orbit_incline*deg_to_rad );
         c_incline = cos( orbit_incline*deg_to_rad );

         x_pos_ECEF(i_sat) = x_pos_orbit*c_long_ascend - y_pos_orbit*c_incline*s_long_ascend;
         y_pos_ECEF(i_sat) = x_pos_orbit*s_long_ascend + y_pos_orbit*c_incline*c_long_ascend; 
         z_pos_ECEF(i_sat) = y_pos_orbit*s_incline;

         los_x_ECEF = x_pos_ECEF(i_sat) - user_x_ECEF; 
         los_y_ECEF = y_pos_ECEF(i_sat) - user_y_ECEF; 
         los_z_ECEF = z_pos_ECEF(i_sat) - user_z_ECEF; 

         rho_user(i_sat) = sqrt ( los_x_ECEF^2 + los_y_ECEF^2 + los_z_ECEF^2 );

         unit_los_x(i_sat) = los_x_ECEF/rho_user(i_sat);
         unit_los_y(i_sat) = los_y_ECEF/rho_user(i_sat);
         unit_los_z(i_sat) = los_z_ECEF/rho_user(i_sat); 
         
         unit_LL_x(i_sat) = CLExx*unit_los_x(i_sat) + CLEyx*unit_los_y(i_sat) + CLEzx*unit_los_z(i_sat);
         unit_LL_y(i_sat) = CLExy*unit_los_x(i_sat) + CLEyy*unit_los_y(i_sat) + CLEzy*unit_los_z(i_sat);
         unit_LL_z(i_sat) = CLExz*unit_los_x(i_sat) + CLEyz*unit_los_y(i_sat) + CLEzz*unit_los_z(i_sat);

         elev_angle(i_sat) = rad_to_deg*asin( up_user_x*unit_los_x(i_sat) + up_user_y*unit_los_y(i_sat) + up_user_z*unit_los_z(i_sat)  );

         i_sat_vis(i_sat) = 0;
         if ( elev_angle(i_sat) > 0.0 ) & ( rho_user(i_sat) < orbit_radius )
            i_sat_vis(i_sat) = 1;
            if ( i_temp < 4 )
               i_temp = i_temp + 1;
% form G matrix from first 4 good satellites for GDOP calculations
               G_matrix(i_temp,1) = unit_LL_x(i_sat);
               G_matrix(i_temp,2) = unit_LL_y(i_sat);
               G_matrix(i_temp,3) = unit_LL_z(i_sat);
               G_matrix(i_temp,4) = 1.0;
               y(i_temp) = elev_angle(i_sat);
               y(i_temp+5) = i_sat;
               y(i_temp+9) = rho_user(i_sat) ;
            end   
         end
        
      end
      
      GDOP_matrix = inv( G_matrix' * G_matrix );
      GDOP = sqrt( GDOP_matrix(1,1) + GDOP_matrix(2,2) + GDOP_matrix(3,3) + GDOP_matrix(4,4) );
      y(5) = GDOP;
      y(14) = t;      
  sys = y;

% end mdlOutputs
