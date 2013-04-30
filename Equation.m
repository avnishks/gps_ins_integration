tiltx = +sin(u(16))*(u(13)*u(20)+u(14)*u(21)+u(15)*u(22))-cos(u(16))*(u(13)*u(17)+u(14)*u(18)+u(15)*u(19));
tiltz = -sin(u(16))*(u( 7)*u(17)+u( 8)*u(18)+u( 9)*u(19))-cos(u(16))*(u( 7)*u(20)+u( 8)*u(21)+u( 9)*u(22));
u(16) = alpha;
% u( 7:15) = Cbn xx xy xz yx yy yz zx zy zz calculated
% u(17:25) = Cbn xx xy xz yx yy yz zx zy zz true

tilty = -(u(10)*u(23)+u(11)*u(24)+u(12)*u(25));



% 
% figure(3)
% llaerror_ins = lla_true - lla_calc;
% poserror_ins = 
% poserror_ins = ecef_true-ecef_calc;
% poserror = [poserror_ins poserror_filt];
% subplot(3,1,1)
% posx     = (poserror(:,1) - poserror(:,7))*.3048;
% cov_posx  = sqrt(filt(:,13)*.3048*.3048);
% plot(t_flight,poserror(:,8),t_flight_samp,cov_posx,'-',t_flight_samp,-cov_posx,'-') %% Later assign variable posx , posy and posz to the first argument
% xlabel('Time (in Seconds)')
% ylabel('Pos Error - x_{ECEF} (m)')
% grid on
% axis auto
% subplot(3,1,2)
% posy     = (poserror(:,2) - poserror(:,8))*.3048;
% cov_posy  = sqrt(filt(:,14)*.3048*.3048);
% plot(t_flight,poserror(:,8),t_flight_samp,cov_posy,'-',t_flight_samp,-cov_posy,'-') %% Later assign variable posx , posy and posz to the first argument
% xlabel('Time (in Seconds)')
% ylabel('Pos Error - y_{ECEF} (m)')
% grid on
% axis auto
% subplot(3,1,3)
% posz     = (poserror(:,3) - poserror(:,9))*.3048;
% cov_posz  = sqrt(filt(:,15)*.3048*.3048);
% plot(t_flight,poserror(:,8),t_flight_samp,cov_posz,'-',t_flight_samp,-cov_posz,'-') %% Later assign variable posx , posy and posz to the first argument
% xlabel('Time (in Seconds)')
% ylabel('Pos Error - z_{ECEF} (m) ')
% grid on
% axis auto