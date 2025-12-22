clear; clc; close all;

%% Parameters of the model
M = 2;      % Mass of the cart [kg]
m1 = 0.2;   % Mass of the first pole [kg]
m2 = 0.2;   % Mass of the second pole [kg]
l1 = 0.6;   % Length of the first pole [m]
l2 = 0.6;   % Length of the second pole [m]
b = 0.1;    % Friction coefficient of the cart [N·s/m]
g = 9.81;   % Gravitational acceleration [m/s^2]

fprintf("=============== Parameters ===============\n")
fprintf("Cart mass: %.2fm\n", M);
fprintf("Pole 1 mass: %.2fm\n", m1);
fprintf("Pole 1 length: %.fm\n", l1);
fprintf("Pole 2 mass: %.2fm\n", m2);
fprintf("Pole 2 length: %.fm\n", l2);
fprintf("Friction coefficient: %.2f\n", b);
fprintf("Gravitational acceleration: %.2fm/s^2\n", g);
fprintf("==========================================\n\n")

%% Linearization

% cosθ₁ ≈ cosθ₂ ≈ 1; cos(θ₁-θ₂) ≈ 1
a11 = M+m1+m2;
a12 = -(m1+m2)*l1;
a13 = -m2*l2;
a21 = a12;
a22 = (m1+m2)*l1^2;
a23 = m2*l1*l2;
a31 = a13;
a32 = a23;
a33 = m2*l2^2;

M0 = [
  a11 a12 a13;
  a21 a22 a23;
  a31 a32 a33
];

invM0 = inv(M0);

B = [1; 0; 0];

D = [b 0 0; zeros(2,3)];

% sin θ₁ ≈ θ₁; sin θ₂ ≈ θ₂; sin(θ₁-θ₂) ≈ θ₁-θ₂
G0 = [
  0 0 0;
  0 -(m1+m2)*l1*g 0;
  0 0 -m2*l2*g
];

%% State-space representation
% State vector: x_state = [x; theta1; theta2; dx; dtheta1; dtheta2]
% Dynamics: qddot = M0^-1 * (B*u - D*qdot - G0*q)
% 
% dx_state/dt = [qdot; qddot] = [qdot; M0^-1*(B*u - D*qdot - G0*q)]
% A * x_state + B_ss * u

A = [
  zeros(3) eye(3);
  -invM0*G0 -invM0*D
];

B_ss = [
  zeros(3,1);
  invM0*B
];

fprintf('State matrix A (6x6):\n');
disp(A);

fprintf('Input matrix B_ss (6x2):\n');
disp(B_ss);

%% Check controllability
fprintf("================ Controllability ================\n")
C_matrix = ctrb(A, B_ss);
rank_C = rank(C_matrix);
fprintf("Controllability matrix rank: %d\n", rank_C);

if rank_C < 6
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

% State weights Q (6x6)
% Order: [x, θ₁, θ₂, dx, dθ₁, dθ₂]
Q = diag([
    1, ...    % x position
    100, ...  % θ₁
    100, ...  % θ₂
    1, ...    % dx
    10, ...   % dθ₁
    0.1       % dθ₂
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

fprintf('\nLQR gain matrix K (1x6):\n');
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
C_pos = [1 0 0 0 0 0];  % Select cart position
k_g = -1 / (C_pos * pinv(Am) * B_ss);

fprintf('Feedforward gain k_g: %.4f\n\n', k_g);

%% Pack parameters for simulation
params = struct("M", M, "m1", m1, ...
  "m2", m2, ...
  "l1", l1, ...
  "l2", l2, ...
  "b", b, ...
  "g", g, ...
  "K_lqr", K_lqr);

%% Initial conditions
x0 = 0.1;             % Initial forward position [m]
theta10 = 0.1;      % Initial angle for pole 1 [rad]
theta20 = -0.3;     % Initial angle for pole 2 [rad]
dx0 = 0;            % Initial forward velocity [m/s]
dtheta10 = 0;       % Initial angular velocity for pole 1 [rad/s]
dtheta20 = 0;       % Initial angular velocity for pole 2 [rad/s]

initial_state = [x0; theta10; theta20; dx0; dtheta10; dtheta20];

fprintf("================= Initial Conditions ==================\n")
fprintf(" Initial cart position: %.2f m\n", x0);
fprintf(" Initial pole 1 angle: %.2f deg\n", rad2deg(theta10));
fprintf(" Initial pole 2 angle: %.2f deg\n", rad2deg(theta20));
fprintf(" Initial cart velocity: %.2f m/s\n", dx0);
fprintf(" Initial pole 1 angular velocity: %.2f deg/s\n", rad2deg(dtheta10));
fprintf(" Initial pole 2 angular velocity: %.2f deg/s\n", rad2deg(dtheta20));
fprintf("=======================================================\n\n")

%% Simulation parameters
sim_time = 100;
tspan = [0 sim_time];

%% Reference signal options:

% Step function
r_step_time = 3.0;  % Step occurs at t=3s
r_amplitude = 0.5;  % Step to 0.5m
r = @(t) (t >= r_step_time) * r_amplitude;

%% LQR control law with reference tracking
% u = -K * x + k_g * r(t)
control_func_lqr = @(t, state) -params.K_lqr * state + k_g * r(t);

%% Solve ODE with LQR control
fprintf("Running control simulation\n\n")
options = odeset("RelTol", 1e-6, "AbsTol", 1e-8);

% ODE dynamics function
f = @(t, y) double_cart_pole_dyn(t, y, params, control_func_lqr);

[t, state] = ode45(f, tspan, initial_state, options);

fprintf("Simulation completed successfully!\\n\\n")

%% Extract states
x = state(:, 1);
theta1 = state(:, 2);
theta2 = state(:, 3);
dx = state(:, 4);
dtheta1 = state(:, 5);
dtheta2 = state(:, 6);

% Compute control inputs over time
u = zeros(length(t), 1);
for i = 1:length(t)
  u(i) = control_func_lqr(t(i), state(i,:)');
end

% Compute state norm (distance from equilibrium)
state_norm = vecnorm(state, 2, 2);

%% Plot results
figure('Position', [100, 100, 1400, 900], 'Color', 'w');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

subplot(3, 3, 1);
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

subplot(3, 3, 2);
plot(t, rad2deg(theta1), 'LineWidth', 2);
hold on;
yline(0, 'k:', 'LineWidth', 1);
hold off;
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\theta_1$ [$^\circ$]', 'FontSize', 14);
title('Pole 1 Angle', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 3);
plot(t, rad2deg(theta2), 'LineWidth', 2);
hold on;
yline(0, 'k:', 'LineWidth', 1);
hold off;
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\theta_2$ [$^\circ$]', 'FontSize', 14);
title('Pole 2 Angle', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 4);
plot(t, dx, 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\dot{x}$ [m/s]', 'FontSize', 14);
title('Cart Velocity', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 5);
plot(t, rad2deg(dtheta1), 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\dot{\theta}_1$ [$^\circ$/s]', 'FontSize', 14);
title('Pole 1 Angular Velocity', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 6);
plot(t, rad2deg(dtheta2), 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$\dot{\theta}_2$ [$^\circ$/s]', 'FontSize', 14);
title('Pole 2 Angular Velocity', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 7);
plot(t, u, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('Force [N]', 'FontSize', 14);
title('Control Input', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 8);
plot(t, rad2deg(theta1), 'b-', 'LineWidth', 1.5, 'DisplayName', '$\theta_1$');
hold on;
plot(t, rad2deg(theta2), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\theta_2$');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
hold off;
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('Angle [$^\circ$]', 'FontSize', 14);
title('Both Pole Angles', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
box off;

subplot(3, 3, 9);
% State norm (distance from equilibrium)
semilogy(t, state_norm, 'b-', 'LineWidth', 2);
xlabel('$t$ [s]', 'FontSize', 14);
ylabel('$||x||_2$', 'FontSize', 14);
title('Distance from Equilibrium', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box off;

sgtitle('LQR-Controlled Double Inverted Pendulum on Cart', 'FontSize', 18, 'FontWeight', 'bold');
