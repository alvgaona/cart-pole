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

invM0 = pinv(M0);

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

fprintf('Input matrix B_ss (6x1):\n');
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

%% Design Reference Model
% State weights Q (6x6)
% Order: [x, θ₁, θ₂, dx, dθ₁, dθ₂]
Q = diag([
    100, ...    % x position
    100, ...  % θ₁
    100, ...  % θ₂
    1, ...    % dx
    10, ...   % dθ₁
    10        % dθ₂ (increased from 0.1 to match dθ₁)
]);

% Input effort weights R (1x1)
R = diag([0.1]);  % Reduced from 1.0 to match single cart-pole

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

Am = A - B_ss * K_lqr;
Bm = B_ss;
Bum = [zeros(4,2); eye(2)];  % Unmatched uncertainty in both pole equations (rows 5,6)

%% Precompute adaptive law matrices
Ts = 0.0001;
phi = pinv(Am) * (expm(Am*Ts) - eye(6)) * [Bm, Bum];
gamma = -pinv(phi) * expm(Am*Ts);

%% Feedforward gain for reference tracking
C_pos = [1 0 0 0 0 0];
k_g = -1/(C_pos * pinv(Am') * Bm);

% Store parameters in a struct
params = struct( ...
  'M', M, ...
  'm1', m1, ...
  'm2', m2, ...
  'l1', l1, ...
  'l2', l2, ...
  'b', b, ...
  'g', g, ...
  'Am', Am, ...
  'Bm', Bm, ...
  'Bum', Bum, ...
  'wf', 20, ...
  'Ts', Ts, ...
  'gamma', gamma, ...
  'k_g', k_g ...
);

params_true = params;
params_true.m1_nominal = m1;
params_true.m1_step_time = 2.0;
params_true.m1_step_value = 0.3;  % m1 increases 50% (multiplicative factor)

%% Simulation parameters
initial_state = [0; 0.1; -0.1; 0; 0; 0];
sim_time = 50;  % Reduced from 100 for easier comparison

[t, state, s_hat, sigma] = cart_pole_l1_solver(...
  @double_cart_pole_dyn, ...
  params_true, ...
  initial_state, ...
  sim_time);

%% Plot results
x = state(:,1); theta1 = state(:,2); theta2 = state(:,3);
dx = state(:,4); dtheta1 = state(:,5); dtheta2 = state(:,6);
x_hat = s_hat(:,1); theta1_hat = s_hat(:,2); theta2_hat = s_hat(:,3);
sigma_m = sigma(:,1); sigma_um1 = sigma(:,2); sigma_um2 = sigma(:,3);

figure('Position', [100, 100, 1500, 900], 'Color', 'w');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

% Cart position
subplot(3,4,1); plot(t, x, 'b', t, x_hat, 'r--', 'LineWidth', 1.5);
legend('Plant', 'Predictor'); title('Cart Position');
xlabel('t [s]'); ylabel('x [m]'); grid on;

% Pole 1 angle
subplot(3,4,2); plot(t, rad2deg(theta1), 'b', t, rad2deg(theta1_hat), 'r--', 'LineWidth', 1.5);
legend('Plant', 'Predictor'); title('Pole 1 Angle');
xlabel('t [s]'); ylabel('$\theta_1$ [deg]'); grid on;

% Pole 2 angle
subplot(3,4,3); plot(t, rad2deg(theta2), 'b', t, rad2deg(theta2_hat), 'r--', 'LineWidth', 1.5);
legend('Plant', 'Predictor'); title('Pole 2 Angle');
xlabel('t [s]'); ylabel('$\theta_2$ [deg]'); grid on;

% Prediction error
subplot(3,4,4); plot(t, vecnorm(state - s_hat, 2, 2), 'LineWidth', 1.5);
title('$||s - \hat{s}||_2$'); xlabel('t [s]'); grid on;

% Cart velocity
subplot(3,4,5); plot(t, dx, 'LineWidth', 1.5);
title('Cart Velocity'); xlabel('t [s]'); ylabel('$\dot{x}$ [m/s]'); grid on;

% Pole 1 angular velocity
subplot(3,4,6); plot(t, rad2deg(dtheta1), 'LineWidth', 1.5);
title('Pole 1 Angular Velocity'); xlabel('t [s]'); ylabel('$\dot{\theta}_1$ [deg/s]'); grid on;

% Pole 2 angular velocity
subplot(3,4,7); plot(t, rad2deg(dtheta2), 'LineWidth', 1.5);
title('Pole 2 Angular Velocity'); xlabel('t [s]'); ylabel('$\dot{\theta}_2$ [deg/s]'); grid on;

% Control input
subplot(3,4,8); plot(t, sigma_m, 'r', 'LineWidth', 2);
title('Control (filtered $-\sigma_m$)'); xlabel('t [s]'); ylabel('Force [N]'); grid on;

% Matched estimate
subplot(3,4,9); plot(t, sigma_m, 'LineWidth', 2);
title('Matched Estimate $\sigma_m$'); xlabel('t [s]'); grid on;

% Unmatched estimate 1
subplot(3,4,10); plot(t, sigma_um1, 'LineWidth', 2);
title('Unmatched Estimate $\sigma_{um,1}$'); xlabel('t [s]'); grid on;

% Unmatched estimate 2
subplot(3,4,11); plot(t, sigma_um2, 'LineWidth', 2);
title('Unmatched Estimate $\sigma_{um,2}$'); xlabel('t [s]'); grid on;

% Reference model poles
subplot(3,4,12);
eigs = eig(params.Am);
plot(real(eigs), imag(eigs), 'x', 'MarkerSize', 10, 'LineWidth', 2);
title('Reference Model Poles'); xlabel('Real'); ylabel('Imag'); grid on;
xlim([-5 5])

sgtitle('Pure L1 Adaptive Control - Double Cart-Pole', 'FontSize', 16);

%% Performance Metrics
fprintf('\n============= Performance Metrics =============\n');

% Position tracking error
x_error = x - 0;  % Reference is 0
rms_x = rms(x_error);
peak_x = max(abs(x_error));
final_x = abs(x_error(end));

% Angle tracking errors
theta1_error = theta1 - 0;
rms_theta1 = rms(rad2deg(theta1_error));
peak_theta1 = max(abs(rad2deg(theta1_error)));
final_theta1 = abs(rad2deg(theta1_error(end)));

theta2_error = theta2 - 0;
rms_theta2 = rms(rad2deg(theta2_error));
peak_theta2 = max(abs(rad2deg(theta2_error)));
final_theta2 = abs(rad2deg(theta2_error(end)));

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

fprintf('\nAngle 1 (theta1):\n');
fprintf('  RMS error:      %.2f deg\n', rms_theta1);
fprintf('  Peak error:     %.2f deg\n', peak_theta1);
fprintf('  Final error:    %.2f deg\n', final_theta1);

fprintf('\nAngle 2 (theta2):\n');
fprintf('  RMS error:      %.2f deg\n', rms_theta2);
fprintf('  Peak error:     %.2f deg\n', peak_theta2);
fprintf('  Final error:    %.2f deg\n', final_theta2);

fprintf('\nPrediction Error:\n');
fprintf('  RMS ||s-s_hat||: %.4f\n', rms_pred);
fprintf('  Peak ||s-s_hat||: %.4f\n', peak_pred);
fprintf('  Final ||s-s_hat||: %.4f\n', final_pred);

fprintf('\nControl Effort:\n');
fprintf('  RMS force:      %.2f N\n', rms_u);
fprintf('  Peak force:     %.2f N\n', peak_u);

fprintf('===============================================\n');
