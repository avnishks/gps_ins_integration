clc
clf
t_flight_samp = linspace(0,600,601);

figure(1)
subplot(3,1,1)
dvx     = (velerror(:,1) - velerror(:,5))*.3048;
cov_vx  = sqrt(filt(:,16)*0.3048*0.3048);
plot(t_flight,dvx,t_flight_samp,cov_vx,'k-',t_flight_samp,-cov_vx,'k-','LineWidth',1)
xlabel('Time (in Seconds)')
ylabel('Velocity Error - V_x (m/s)')
grid on
h = legend('Estimation Error','$\sqrt{P_{44}}$');
set(h,'Interpreter','latex')
axis auto
subplot(3,1,2)
dvy     = (velerror(:,2) - velerror(:,4))*.3048;
cov_vy  = sqrt(filt(:,17)*0.3048*0.3048);
plot(t_flight,dvy,t_flight_samp,cov_vy,'k-',t_flight_samp,-cov_vy,'k-')
xlabel('Time (in Seconds)')
ylabel('Velocity Error - V_y (m/s)')
grid on
h = legend('Estimation Error','$\sqrt{P_{55}}$');
set(h,'Interpreter','latex')
axis auto
subplot(3,1,3)
dvz     = (velerror(:,3) - velerror(:,6))*.3048;
cov_vz  = sqrt(filt(:,18)*0.3048*0.3048);
plot(t_flight,dvz',t_flight_samp,cov_vz,'k-',t_flight_samp,-cov_vz,'k-')
xlabel('Time (in Seconds)')
ylabel('Velocity Error - V_z (m/s')
grid on
axis auto
h = legend('Estimation Error','$\sqrt{P_{66}}$');
set(h,'Interpreter','latex')

% figure(101)
% subplot(2,1,1)
% vxtrue = truevel(:,1)*0.3048;
% vxins  = navvel(:,4)*0.3048;
% vxfilt = vxins + velerror(:,5)*0.3048;
% plot(t_flight,vxtrue,t_flight,vxins,t_flight,vxfilt,'LineWidth',1)
% grid on
% xlabel('Time (in Seconds)')
% ylabel('Velocity - V_x (m/s)')
% legend('True','INS','Filtered')
% subplot(2,1,2)
% plot(t_flight,vxtrue - vxins,t_flight,vxtrue - vxfilt,'LineWidth',1)
% grid on
% xlabel('Time (in Seconds)')
% ylabel('Velocity Error- V_x (m/s)')
% legend('INS','Filtered')

%%
figure(2)
subplot(3,1,1)
tiltx     = (tilterror(:,1) - tilterror(:,4))*1;
cov_tiltx  = sqrt(filt(:,19)*1*1);
plot(t_flight,tiltx,t_flight_samp,cov_tiltx,'k-',t_flight_samp,-cov_tiltx,'k-')
xlabel('Time (in Seconds)')
ylabel('Tilt Error - x (in rad)')
grid on
h = legend('Estimation Error','$\sqrt{P_{77}}$');
set(h,'Interpreter','latex')
axis auto
subplot(3,1,2)
tilty     = (tilterror(:,2) - tilterror(:,5))*1;
cov_tilty  = sqrt(filt(:,20)*1*1);
plot(t_flight,tilty,t_flight_samp,cov_tilty,'k-',t_flight_samp,-cov_tilty,'k-')
xlabel('Time (in Seconds)')
ylabel('Tilt Error - y (in rad)' )
grid on
h = legend('Estimation Error','$\sqrt{P_{88}}$');
set(h,'Interpreter','latex')
axis auto
subplot(3,1,3)
tiltz     = (tilterror(:,3) - tilterror(:,6))*1;
cov_tiltz  = sqrt(filt(:,21)*1*1);
plot(t_flight,tiltz,t_flight_samp,cov_tiltz,'k-',t_flight_samp,-cov_tiltz,'k-')
xlabel('Time (in Seconds)')
ylabel('Tilt Error - z (in rad)')
grid on
h = legend('Estimation Error','$\sqrt{P_{99}}$');
set(h,'Interpreter','latex')
axis auto
%%
figure(3)
subplot(3,1,1)
plot(lla_true.signals.values(:,2)/0.01745329,lla_true.signals.values(:,1)/0.01745329,lla_calc.signals.values(:,2),lla_calc.signals.values(:,1),'LineWidth',1)
grid on
xlabel('Longitude (in degrees)')
ylabel('Latitude (in degrees)')
axis auto
subplot(3,1,2)
plot(lla_true.time,lla_true.signals.values(:,1)/0.01745329-lla_calc.signals.values(:,1),'LineWidth',2)
grid on
xlabel('Time (in Seconds)')
ylabel('Error in Latitude (in degrees)')
axis auto
subplot(3,1,3)
plot(lla_true.time,lla_true.signals.values(:,2)/0.01745329-lla_calc.signals.values(:,2),'LineWidth',2)
grid on
xlabel('Time (in Seconds)')
ylabel('Error in Longitude (in degrees)')
axis auto

figure(103)
plot(lla_true.signals.values(:,2)/0.01745329,lla_true.signals.values(:,1)/0.01745329,lla_calc.signals.values(:,2),lla_calc.signals.values(:,1),'LineWidth',2)
grid on
xlabel('Longitude (in degrees)')
ylabel('Latitude (in degrees)')
legend('True Trajectory',' INS ')
axis equal

%% 
% Other Plots Required:
%     VELOCITY (NOT VELOCITY ERROR)
%     POSITION ERROR ESTIMATE
%     ANGULAR RATE (BODY)
%     GYRO MEASUREMENT
%     ACCELEROMETER MEASUREMENT
%     NOISE OF BOTH GYRO AND ACCELEROMETER
%     PSEUDORANGE MEASUREMENT
%     PSEUDORANGE NOISE
%     CLOCK BIAS
%     CLOCK DRIFT
%     
    
%%
figure(4)
poserror_ins = ecef_true.signals.values-ecef_calc.signals.values;
poserror = [poserror_ins poserror_filt];
cons = 0.3048;
ka = 20;
subplot(3,1,1)
posx     = (poserror_filt(:,1) + poserror_filt(:,4))*.3048;
cov_posx  = sqrt(filt(:,13)*.3048*.3048);
plot(t_flight,posx*cons,t_flight_samp,cov_posx,'k-',t_flight_samp,-cov_posx,'k-') %% Later assign variable posx , posy and posz to the first argument
xlabel('Time (in Seconds)')
ylabel('Pos Error - x_{ECEF} (m)')
grid on
h = legend('Estimation Error','$\sqrt{P_{11}}$');
set(h,'Interpreter','latex')
axis([0 600 -ka ka])
subplot(3,1,2)
posy     = (poserror_filt(:,2) + poserror_filt(:,5))*.3048;
cov_posy  = sqrt(filt(:,14)*.3048*.3048);
plot(t_flight,posy*cons,t_flight_samp,cov_posy,'k-',t_flight_samp,-cov_posy,'k-') %% Later assign variable posx , posy and posz to the first argument
xlabel('Time (in Seconds)')
ylabel('Pos Error - y_{ECEF} (m)')
h = legend('Estimation Error','$\sqrt{P_{22}}$');
set(h,'Interpreter','latex')
grid on
axis([0 600 -ka ka])
subplot(3,1,3)
posz     = (poserror_filt(:,3) + poserror_filt(:,6))*.3048;
cov_posz  = sqrt(filt(:,15)*.3048*.3048);
plot(t_flight,posz*cons,t_flight_samp,cov_posz,'k-',t_flight_samp,-cov_posz,'k-') %% Later assign variable posx , posy and posz to the first argument
xlabel('Time (in Seconds)')
ylabel('Pos Error - z_{ECEF} (m) ')
grid on
h = legend('Estimation Error','$\sqrt{P_{33}}$');
set(h,'Interpreter','latex')
axis([0 600 -ka ka])
%% 
% figure(999)
% subplot(3,1,1)
% plot(t_flight,abs(posx4)*0.3048,t_flight,abs(posx3(1:18002,1))*0.3048,t_flight,abs(posx1(1:18002,1))*0.3048)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Position Error (Absolute Value) - x (m)')
% legend('4 Satellites','3 Satellites','1 Satellite')
% title('Position Estimates - Comparison of Visible Satellites')
% subplot(3,1,2)
% plot(t_flight,abs(posy4)*0.3048,t_flight,abs(posy3(1:18002,1))*0.3048,t_flight,abs(posy1(1:18002,1))*0.3048)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Position Error(Absolute Value) - y (m)')
% legend('4 Satellites','3 Satellites','1 Satellite')
% subplot(3,1,3)
% plot(t_flight,abs(posz4)*0.3048,t_flight,abs(posz3(1:18002,1))*0.3048,t_flight,abs(posz1(1:18002,1))*0.3048)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Position Error (Absolute Value) - z (m)')
% legend('4 Satellites','3 Satellites','1 Satellite')
% 
% figure(997)
% subplot(3,1,1)
% plot(t_flight,abs(dvx4)*1,t_flight,abs(dvx3(1:18002,1))*1,t_flight,abs(dvx1(1:18002,1))*1)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Velocity Error (Absolute Value) - x (m/s)')
% legend('4 Satellites','3 Satellites','1 Satellite')
% title('Absolute ValueVelocity Estimates - Comparison of Visible Satellites')
% subplot(3,1,2)
% plot(t_flight,abs(dvy4)*1,t_flight,abs(dvy3(1:18002,1))*1,t_flight,abs(dvy1(1:18002,1))*1)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Velocity Error(Absolute Value) - y (m/s)')
% legend('4 Satellites','3 Satellites','1 Satellite')
% subplot(3,1,3)
% plot(t_flight,abs(dvz4)*1,t_flight,abs(dvz3(1:18002,1))*1,t_flight,abs(dvz1(1:18002,1))*1)
% grid on
% xlabel('Time (in seconds)')
% ylabel('Velocity Error (Absolute Value) - z (m/s)')
% legend('4 Satellites','3 Satellites','1 Satellite')
%%
figure(5)
bank = [0,0,25.,25.,0,0,-25.,-25.,0,0,0,0,0,0,0,0,0,0];
tii = [0,60.,70.,90.,100.,160.,170.,190.,200.,260.,270.,290.,300.,360.,370.,390.,400.,600.];
hold on
plot(tii,bank,'LineWidth',2)
grid on
xlabel('Time (in Seconds)')
ylabel('Aircraft Bank Angle (in degrees)')
title('Aircraft Bank Angle')
axis([0 600 -30 30])
%%
figure(6)
plot(t_flight_samp,rho_noise(1:601,4),'LineWidth',1.5)
grid on
xlabel('Time (in Seconds)')
ylabel('Noise In Pseudorange Measurement (m)')
title('Pseudorange Measurement')
axis auto
%%
figure(7)
subplot(2,1,1)
plot(t_flight,ideal_meas(:,1),t_flight,noisy_meas(:,1))
grid on
xlabel('Time (in Seconds)')
ylabel('\omega_1 (in rad/sec)')
legend('True','Measurement')
title('Gyroscope Measurement ')
subplot(2,1,2)
plot(t_flight,(ideal_meas(:,1)-noisy_meas(:,1))*(180/pi))
grid on
xlabel('Time (in Seconds)')
ylabel('Error in \omega_1 (in degrees/sec)')
%%
figure(8)
subplot(2,1,1)
plot(t_flight,ideal_meas(:,6)*0.3048,t_flight,noisy_meas(:,6)*0.3048)
grid on
xlabel('Time (in Seconds)')
ylabel('a_3(in m/sec^2)')
legend('True','Measurement')
title('Accelerometer Measurement ')
subplot(2,1,2)
plot(t_flight,(ideal_meas(:,4)-noisy_meas(:,4))*0.3048)
grid on
xlabel('Time (in Seconds)')
ylabel('Error in a_3 (in m/sec^2)')
%%
% figure(9)
% plot(t_flight,filt1,'LineWidth',2)
% grid on
% xlabel('Time (in Seconds)')
% ylabel('GDOP')