clear; clc; close all;

%% Parameters of the model
M = 2;      % Mass of the cart [kg]
m = 0.2;    % Mass of the pole [kg]
l = 0.6;    % Length of the first pole [m]
b = 0.1;    % Friction coefficient of the cart [N·s/m]
g = 9.81;   % Gravitational acceleration [m/s^2]

fprintf("=============== Parameters ===============\n")
fprintf("Cart mass: %.2fm\n", M);
fprintf("Pole mass: %.2fm\n", m);
fprintf("Pole length: %.fm\n", l);
fprintf("Friction coefficient: %.2f\n", b);
fprintf("Gravitational acceleration: %.2fm/s^2\n", g);
fprintf("==========================================\n\n")

%% Linearization

% cosθ ≈ 1
a11 = M+m;
a12 = -m*l;
a21 = a12;
a22 = m*l^2;

M0 = [
  a11 a12;
  a21 a22
];

invM0 = inv(M0);

B = [1; 0];

D = [b 0; 0 0];

% sin θ ≈ θ
G0 = [
  0 0;
  0 -m*l*g;
];

%% State-space representation
% State vector: x_state = [x; theta; dx; dtheta]
% Dynamics: qddot = M0^-1 * (B*u - D*qdot - G0*q)
% 
% dx_state/dt = [qdot; qddot] = [qdot; M0^-1*(B*u - D*qdot - G0*q)]
% A * x_state + B_ss * u

A = [
  zeros(2) eye(2);
  -invM0*G0 -invM0*D
];

B_ss = [
  zeros(2,1);
  invM0*B
];

fprintf('State matrix A (4x4):\n');
disp(A);

fprintf('Input matrix B_ss (4x1):\n');
disp(B_ss);

%% Check controllability
fprintf("================ Controllability ================\n")
C_matrix = ctrb(A, B_ss);
rank_C = rank(C_matrix);
fprintf("Controllability matrix rank: %d\n", rank_C);

if rank_C < 4
  warning(" System is not fully controllable!")
else
  fprintf(" System is fully controllable.\n")
end
fprintf("=================================================\n\n")

%% Check open-loop eigenvalues
fprintf("============ Open-loop Eigenvalues ==============\n")
open_eigval = eig(A);
fprintf("Open-loop eigenvalues:\n")
disp(open_eigval)

unstable_modes = sum(real(open_eigval) > 0);
fprintf("Number of unstable modes: %d\n", unstable_modes);
fprintf("=================================================\n")

%% Design LQR Controller

% State weights Q (4x4)
% Order: [x, θ, dx, dθ]
Q = diag([
    1, ...    % x position
    100, ...  % θ
    1, ...    % dx
    10, ...   % dθ
]);

% Input effort weights R (1x1)
R = diag([0.1]);

fprintf('\n===== Designing LQR Controller ====\n');
fprintf('State weights Q (diagonal):\n');
disp(diag(Q)');
fprintf('Control weights R (diagonal):\n');
disp(diag(R)');

% Compute LQR gain
[K_lqr, S, poles_closed] = lqr(A, B_ss, Q, R);

fprintf('\nLQR gain matrix K (1x4):\n');
disp(K_lqr);

fprintf('\nClosed-loop poles:\n');
disp(poles_closed);

if all(real(poles_closed) < 0)
  fprintf("All closed-loop poles have negative real parts.\n\n");
else
  fprintf("Some closed-loop poles have positive real part.\n\n")
end

%% Feedforward gain for reference tracking
% Compute DC gain from control input to cart position (first state)
Am = A - B_ss * K_lqr;
C_pos = [1 0 0 0];  % Select cart position
k_g = -1 / (C_pos * pinv(Am) * B_ss);

fprintf('Feedforward gain k_g: %.4f\n\n', k_g);

%% Pack parameters for simulation
params = struct("M", M, "m", m, ...
  "l", l, ...
  "b", b, ...
  "g", g, ...
  "K_lqr", K_lqr);

%% Initial conditions
x0 = 0;             % Initial forward position [m]
theta0 = 0.1;       % Initial angle for pole [rad]
dx0 = 0;            % Initial forward velocity [m/s]
dtheta0 = 0;        % Initial angular velocity for pole [rad/s]

initial_state = [x0; theta0; dx0; dtheta0];

fprintf("================= Initial Conditions ==================\n")
fprintf(" Initial cart position: %.2f m\n", x0);
fprintf(" Initial pole angle: %.2f deg\n", rad2deg(theta0));
fprintf(" Initial cart velocity: %.2f m/s\n", dx0);
fprintf(" Initial pole angular velocity: %.2f deg/s\n", rad2deg(dtheta0));
fprintf("=======================================================\n\n")

%% Simulation parameters
sim_time = 30;
tspan = [0 sim_time];

%% Reference signal: Step function
r_step_time = 3.0;  % Step occurs at t=3s
r_amplitude = 0.5;  % Step to 0.5m
r = @(t) (t >= r_step_time) * r_amplitude;

%% LQR control law with feedforward
% u = -K_lqr * x + k_g * r(t)

control_func_lqr = @(t, state) -params.K_lqr * state + k_g * r(t);

%% Solve ODE with LQR control
fprintf("Running LQR control simulation\n\n")
options = odeset("RelTol", 1e-6, "AbsTol", 1e-8);

% ODE dynamics function
f = @(t, y) single_cart_pole_dyn(t, y, params, control_func_lqr);

[t, state] = ode45(f, tspan, initial_state, options);

fprintf("Simulation completed successfully!\n\n")

%% Extract states
x = state(:, 1);
theta = state(:, 2);
dx = state(:, 3);
dtheta = state(:, 4);

% Compute control inputs over time
u = zeros(length(t), 1);
for i = 1:length(t)
  u(i) = control_func_lqr(t(i), state(i,:)');
end

% Compute state norm (distance from equilibrium)
state_norm = vecnorm(state, 2, 2);

%% Plot results
figure('Position', [100, 100, 1200, 800], 'Color', 'w');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

subplot(2, 3, 1);
r_vec = arrayfun(r, t);
plot(t, x, 'b-', 'LineWidth', 2, 'DisplayName', 'Plant');
hold on;
plot(t, r_vec, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold off;
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$x$ [m]', 'FontSize', 14);
title('Cart Position', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(2, 3, 2);
plot(t, rad2deg(theta), 'LineWidth', 2);
hold on;
yline(0, 'k:', 'LineWidth', 1);
hold off;
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\theta$ [$^\circ$]', 'FontSize', 14);
title('Pole Angle', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(2, 3, 3);
plot(t, dx, 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14);
title('Cart Velocity', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(2, 3, 4);
plot(t, rad2deg(dtheta), 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\dot{\theta}$ [$^\circ$/s]', 'FontSize', 14);
title('Pole Angular Velocity', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(2, 3, 5);
plot(t, u, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('Force [N]', 'FontSize', 14);
title('Control Input', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(2, 3, 6);
% State norm (distance from equilibrium)
semilogy(t, state_norm, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$||x||_2$', 'FontSize', 14);
title('Distance from Equilibrium', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

sgtitle('LQR-Controlled Single Inverted Pendulum on Cart', 'FontSize', 18, 'FontWeight', 'bold');
