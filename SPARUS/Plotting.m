%% here you can implement the code in order to have some figures ...
%%
% x_max = max(PosE_S(:,1));
% y_max = max(PosE_S(:,2));
% z_max = max(PosE_S(:,3));
% phi_max = max(PosE_S(:,4));
% theta_max = max(PosE_S(:,5));
% psi_max = max(PosE_S(:,6));
% 
% ylim_linear = max([x_max y_max z_max])
% ylim_angular = max([phi_max theta_max psi_max])
%% Position
figure('Name','Simulated Position');
subplot(6,1,1);
plot(PosE_S(:,1))
title('x Comparison')
xlabel('time (s)')
ylabel('x (m)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,2);
plot(PosE_S(:,2))
title('y Comparison')
xlabel('time (s)')
ylabel('y (m)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,3);
plot(PosE_S(:,3))
title('z Comparison')
xlabel('time (s)')
ylabel('z (m)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,4);
plot(PosE_S(:,4))
title('phi Comparison')
xlabel('time (s)')
ylabel('phi (rad)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,5);
plot(PosE_S(:,5))
title('theta Comparison')
xlabel('time (s)')
ylabel('theta (rad)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,6);
plot(PosE_S(:,6))
title('psi Comparison')
xlabel('time (s)')
ylabel('psi (rad)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

%% Velocity
figure('Name','Simulated Velocity');
subplot(6,1,1);
plot(VitB_S(:,1))
title('u Comparison')
xlabel('time (s)')
ylabel('u (m/s)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,2);
plot(VitB_S(:,2))
title('v Comparison')
xlabel('time (s)')
ylabel('v (m/s)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,3);
plot(VitB_S(:,3))
title('w Comparison')
xlabel('time (s)')
ylabel('w (m/s)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,4);
plot(VitB_S(:,4))
title('p Comparison')
xlabel('Time')
ylabel('p (rad/s)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,5);
plot(VitB_S(:,5))
title('q Comparison')
xlabel('time (s)')
ylabel('q (rad/s)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,6);
plot(VitB_S(:,6))
title('r Comparison')
xlabel('time (s)')
ylabel('r (rad/s)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

%% Acceleration
figure('Name','Simulated Acceleration');
subplot(6,1,1);
plot(AccB_S(:,1))
title('udot Comparison')
xlabel('time (s)')
ylabel('udot (m/s2)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,2);
plot(AccB_S(:,2))
title('vdot Comparison')
xlabel('time (s)')
ylabel('vdot (m/s2)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,3);
plot(AccB_S(:,3))
title('wdot Comparison')
xlabel('time (s)')
ylabel('wdot (m/s2)')
% ylim([-ylim_linear ylim_linear])
hold on
grid on

subplot(6,1,4);
plot(AccB_S(:,4))
title('pdot Comparison')
xlabel('Time')
ylabel('p (rad/s2)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,5);
plot(AccB_S(:,5))
title('qdot Comparison')
xlabel('time (s)')
ylabel('q (rad/s2)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

subplot(6,1,6);
plot(AccB_S(:,6))
title('rdot Comparison')
xlabel('time (s)')
ylabel('r (rad/s2)')
% ylim([-ylim_angular ylim_angular])
hold on
grid on

%% Save
filename1 = 'WithAppendages.mat';
filename2 = 'WithoutAppendages.mat';

save(filename2, 'PosE_S', 'VitB_S', 'AccB_S', 'Thrust_S')


