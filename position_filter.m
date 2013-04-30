deg_to_rad = 0.01745329252;
ellipticity = 3.35281066474E-03;
R_equator = 20925604.0;

for i = 1:length(lla_true)
    user_lat  = lla_true(i,1);
    user_long = lla_true(i,2);
    user_alt  = lla_true(i,3);
    s_user_lat = sin ( user_lat*1 );
    c_user_lat = cos ( user_lat*1 );
    s_user_long = sin ( user_long*1 );
    c_user_long = cos ( user_long*1 );
    sqr_one_f = ( 1.0 - ellipticity )*( 1.0 - ellipticity );
    denom = sqrt ( (c_user_lat*c_user_lat) + sqr_one_f*(s_user_lat*s_user_lat) );
    R = R_equator/denom + user_alt;
    S = R_equator*sqr_one_f/denom + user_alt;
    user_x_ECEF(i) = R*c_user_lat*c_user_long;
    user_y_ECEF(i) = R*c_user_lat*s_user_long;
    user_z_ECEF(i) = S*s_user_lat;
end
true_ECEF_man = [user_x_ECEF; user_y_ECEF; user_z_ECEF]';

for i = 1:length(lla_calc)
    user_lat  = lla_calc(i,1);
    user_long = lla_calc(i,2);
    user_alt  = lla_calc(i,3);
    s_user_lat = sin ( user_lat*deg_to_rad );
    c_user_lat = cos ( user_lat*deg_to_rad );
    s_user_long = sin ( user_long*deg_to_rad );
    c_user_long = cos ( user_long*deg_to_rad );
    sqr_one_f = ( 1.0 - ellipticity )*( 1.0 - ellipticity );
    denom = sqrt ( (c_user_lat*c_user_lat) + sqr_one_f*(s_user_lat*s_user_lat) );
    R = R_equator/denom + user_alt;
    S = R_equator*sqr_one_f/denom + user_alt;
    user_x_ECEF(i) = R*c_user_lat*c_user_long;
    user_y_ECEF(i) = R*c_user_lat*s_user_long;
    user_z_ECEF(i) = S*s_user_lat;
    
end
calc_ECEF_man = [user_x_ECEF ;user_y_ECEF ;user_z_ECEF]';