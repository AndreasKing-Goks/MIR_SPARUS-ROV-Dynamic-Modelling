% Load Data
WithoutAppendages = load('WithoutAppendages.mat');
WithAppendages = load('WithAppendages.mat');

%% Position
figure('Name','Simulated Position');
subplot(6,1,1); 
plot(WithAppendages.PosE_S(:,1))
hold on
plot(WithoutAppendages.PosE_S(:,1),'r')
title('x Comparison')
xlabel('time (s)')
ylabel('x (m)')
xlim([0 1000])
legend("With Appendages", "Without Appendages", "Location", "northwest")
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,2);
plot(WithAppendages.PosE_S(:,2))
hold on
plot(WithoutAppendages.PosE_S(:,2),'r')
title('y Comparison')
xlabel('time (s)')
ylabel('y (m)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,3);
plot(WithAppendages.PosE_S(:,3))
hold on
plot(WithoutAppendages.PosE_S(:,3),'r')
title('z Comparison')
xlabel('time (s)')
ylabel('z (m)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,4);
plot(WithAppendages.PosE_S(:,4))
hold on
plot(WithoutAppendages.PosE_S(:,4),'r')
title('phi Comparison')
xlabel('time (s)')
ylabel('phi (rad)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,5);
plot(WithAppendages.PosE_S(:,5))
hold on
plot(WithoutAppendages.PosE_S(:,5),'r')
title('theta Comparison')
xlabel('time (s)')
ylabel('theta (rad)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,6);
plot(WithAppendages.PosE_S(:,6))
hold on
plot(WithoutAppendages.PosE_S(:,6),'r')
title('psi Comparison')
xlabel('time (s)')
ylabel('psi (rad)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

%% Velocity
figure('Name','Simulated Velocity');
subplot(6,1,1);
plot(WithAppendages.VitB_S(:,1))
hold on
plot(WithoutAppendages.VitB_S(:,1),'r')
title('u Comparison')
xlabel('time (s)')
ylabel('u (m/s)')
xlim([0 1000])
legend("With Appendages", "Without Appendages", "Location", "northwest")
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,2);
plot(WithAppendages.VitB_S(:,2))
hold on
plot(WithoutAppendages.VitB_S(:,2),'r')
title('v Comparison')
xlabel('time (s)')
ylabel('v (m/s)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,3);
plot(WithAppendages.VitB_S(:,3))
hold on
plot(WithoutAppendages.VitB_S(:,3),'r')
title('w Comparison')
xlabel('time (s)')
ylabel('w (m/s)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,4);
plot(WithAppendages.VitB_S(:,4))
hold on
plot(WithoutAppendages.VitB_S(:,4),'r')
title('p Comparison')
xlabel('Time')
ylabel('p (rad/s)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,5);
plot(WithAppendages.VitB_S(:,5))
hold on
plot(WithoutAppendages.VitB_S(:,5),'r')
title('q Comparison')
xlabel('time (s)')
ylabel('q (rad/s)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,6);
plot(WithAppendages.VitB_S(:,6))
hold on
plot(WithoutAppendages.VitB_S(:,6),'r')
title('r Comparison')
xlabel('time (s)')
ylabel('r (rad/s)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

%% Acceleration
figure('Name','Simulated Acceleration');
subplot(6,1,1);
plot(WithAppendages.AccB_S(:,1))
hold on
plot(WithoutAppendages.AccB_S(:,1),'r')
title('udot Comparison')
xlabel('time (s)')
ylabel('udot (m/s2)')
xlim([0 1000])
legend("With Appendages", "Without Appendages", "Location", "northwest")
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,2);
plot(WithAppendages.AccB_S(:,2))
hold on
plot(WithoutAppendages.AccB_S(:,2),'r')
title('vdot Comparison')
xlabel('time (s)')
ylabel('vdot (m/s2)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,3);
plot(WithAppendages.AccB_S(:,3))
hold on
plot(WithoutAppendages.AccB_S(:,3),'r')
title('wdot Comparison')
xlabel('time (s)')
ylabel('wdot (m/s2)')
xlim([0 1000])
% ylim([-ylim_linear ylim_linear])
grid on

subplot(6,1,4);
plot(WithAppendages.AccB_S(:,4))
hold on
plot(WithoutAppendages.AccB_S(:,4),'r')
title('pdot Comparison')
xlabel('Time')
ylabel('p (rad/s2)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,5);
plot(WithAppendages.AccB_S(:,5))
hold on
plot(WithoutAppendages.AccB_S(:,5),'r')
title('qdot Comparison')
xlabel('time (s)')
ylabel('q (rad/s2)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on

subplot(6,1,6);
plot(WithAppendages.AccB_S(:,6))
hold on
plot(WithoutAppendages.AccB_S(:,6),'r')
title('rdot Comparison')
xlabel('time (s)')
ylabel('r (rad/s2)')
xlim([0 1000])
% ylim([-ylim_angular ylim_angular])
grid on