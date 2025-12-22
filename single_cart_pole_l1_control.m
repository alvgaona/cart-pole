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

invM0 = pinv(M0);

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

%% Design Reference Model
% State weights Q (4x4)
% Order: [x, θ, dx, dθ]
Q = diag([
    100, ...    % x position
    100, ...  % θ₁
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

Am = A - B_ss * K_lqr;
Bm = B_ss;
Bum = [0; 0; 0; 1];

%% Precompute adaptive law matrices
Ts = 0.0001;
phi = pinv(Am) * (expm(Am*Ts) - eye(4)) * [Bm, Bum];
gamma = -pinv(phi) * expm(Am*Ts);

%% Feedforward gain for reference tracking
% Compute DC gain from control input to cart position (first state)
C_pos = [1 0 0 0];  % Select cart position
k_g = -1 / (C_pos * pinv(Am) * Bm);

% Store parameters in a struct
params = struct( ...
  'M', M, ...
  'm', m, ...
  'l', l, ...
  'b', b, ...
  'g', g, ...
  'Am', Am, ...
  'Bm', Bm, ...
  'Bum', Bum, ...
  'K_lqr', K_lqr, ...
  'wf', 20, ...
  'Ts', Ts, ...
  'gamma', gamma, ...
  'k_g', k_g ...
);

params_true = params;

% Mass uncertainty
params_true.m_nominal = m;
params_true.m_step_time = 2.0;
params_true.m_step_value = 0.3;  % m: 0.2 → 0.3 kg (+50%)

% Friction uncertainty
params_true.b_nominal = b;
params_true.b_step_time = 5.0;
params_true.b_step_value = 0.2;  % b: 0.1 → 0.2 N·s/m (+100%)

% Length uncertainty
params_true.l_nominal = l;
params_true.l_step_time = 8.0;
params_true.l_step_value = 0.7;  % l: 0.6 → 0.7 m (+17%)

%% Simulation parameters
initial_state = [0; 0.6; 0; 0];
sim_time = 50;

% Step function
r_step_time = 3.0;  % Step occurs at t=3s
r_amplitude = 0.5;  % Step to 0.5m
r = @(t) (t >= r_step_time) * r_amplitude;

%% Disturbance options (returns [F_cart; tau_pole])

% No disturbance (default)
% dist = @(t, x) [0; 0];

% Impulse force on cart (5N for 0.1s at t=10s)
% dist = @(t, x) [(t >= 10.0 && t < 10.1) * 5.0; 0];

% Step force on cart (2N starting at t=7s)
% dist = @(t, x) [(t >= 7.0) * 2.0; 0];

% Sinusoidal force on cart (1N at 1Hz)
% dist = @(t, x) [sin(2*pi*1*t); 0];

% Impulse torque on pole (0.5 N·m for 0.1s at t=12s)
% dist = @(t, x) [0; (t >= 12.0 && t < 12.1) * 0.5];

% Combined: step force + sinusoidal torque
dist = @(t, x) [(t >= 7.0) * 2.0; 0.01*sin(2*pi*0.5*t)];

[t, state, s_hat, sigma] = cart_pole_l1_solver(...
  @single_cart_pole_dyn, ...
  params_true, ...
  initial_state, ...
  sim_time, ...
  r, ...
  dist);

%% Plot results
x = state(:,1); theta = state(:,2); dx = state(:,3); dtheta = state(:,4);
x_hat = s_hat(:,1); theta_hat = s_hat(:,2);
sigma_m = sigma(:,1); sigma_um = sigma(:,2);

figure('Position', [100, 100, 1500, 600], 'Color', 'w');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% States
subplot(2,4,1);
r_vec = arrayfun(r, t)';
plot(t, x, 'b', t, x_hat, 'r--', t, r_vec, 'k:', 'LineWidth', 1.5);
legend('Plant', 'Predictor', 'Reference'); title('Cart Position');
xlabel('t [s]'); ylabel('x [m]'); grid on;

subplot(2,4,2); plot(t, rad2deg(theta), 'b', t, rad2deg(theta_hat), 'r--', 'LineWidth', 1.5);
legend('Plant', 'Predictor'); title('Pole Angle');
xlabel('t [s]'); ylabel('$\theta$ [deg]'); grid on;

% Prediction error
subplot(2,4,3); plot(t, abs(state - s_hat), 'LineWidth', 1.5);
title('|s - $\hat{s}$|'); xlabel('t [s]'); grid on;

% Control input
subplot(2,4,4); plot(t, sigma_m, 'r', 'LineWidth', 2);
title('Control (filtered - $\sigma_m$)'); xlabel('t [s]'); ylabel('Force [N]'); grid on;

% Uncertainty estimates
subplot(2,4,5); plot(t, sigma_m, 'LineWidth', 2);
title('Matched Estimate $\sigma_m$'); xlabel('t [s]'); grid on;

subplot(2,4,6); plot(t, sigma_um, 'LineWidth', 2);
title('Unmatched Estimate $\sigma_{um}$'); xlabel('t [s]'); grid on;

% Eigenvalues of reference model
subplot(2,4,7);
eigs = eig(params.Am);
plot(real(eigs), imag(eigs), 'x', 'MarkerSize', 10, 'LineWidth', 2);
title('Reference Model Poles'); xlabel('Real'); ylabel('Imag'); grid on;
xlim([-5 5])

sgtitle('L1 Adaptive Control Results', 'FontSize', 16);

%% Performance Metrics
fprintf('\n============= Performance Metrics =============\n');

% Position tracking error
r_vec = arrayfun(r, t)';  % Evaluate reference over time (column vector)
x_error = x - r_vec;
rms_x = rms(x_error);
peak_x = max(abs(x_error));
final_x = abs(x_error(end));

% Angle tracking errors
theta_error = theta - 0;
rms_theta = rms(rad2deg(theta_error));
peak_theta = max(abs(rad2deg(theta_error)));
final_theta = abs(rad2deg(theta_error(end)));

% Prediction error
pred_error = vecnorm(state - s_hat, 2, 2);
rms_pred = rms(pred_error);
peak_pred = max(pred_error);
final_pred = pred_error(end);

% Control effort
u = sigma_m;
rms_u = rms(u);
peak_u = max(abs(u));

% Settling time (2% criterion for position)
settling_threshold = 0.02 * max(abs(x));
settled_idx = find(abs(x_error) > settling_threshold, 1, 'last');
if isempty(settled_idx)
    settling_time = t(1);
else
    settling_time = t(min(settled_idx + 1, length(t)));
end

fprintf('Position (x):\n');
fprintf('  RMS error:      %.4f m\n', rms_x);
fprintf('  Peak error:     %.4f m\n', peak_x);
fprintf('  Final error:    %.4f m\n', final_x);
fprintf('  Settling time:  %.2f s\n', settling_time);

fprintf('\nAngle (theta):\n');
fprintf('  RMS error:      %.2f deg\n', rms_theta);
fprintf('  Peak error:     %.2f deg\n', peak_theta);
fprintf('  Final error:    %.2f deg\n', final_theta);

fprintf('\nPrediction Error:\n');
fprintf('  RMS ||s-s_hat||: %.4f\n', rms_pred);
fprintf('  Peak ||s-s_hat||: %.4f\n', peak_pred);
fprintf('  Final ||s-s_hat||: %.4f\n', final_pred);

fprintf('\nControl Effort:\n');
fprintf('  RMS force:      %.2f N\n', rms_u);
fprintf('  Peak force:     %.2f N\n', peak_u);

fprintf('===============================================\n');




