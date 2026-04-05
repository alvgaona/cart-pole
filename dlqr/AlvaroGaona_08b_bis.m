clear; clc; close all;

%% Student parameters
numero = 8;
rng(numero);
valor = rand;

M = 1 + valor/2;
m = 0.2 * (1 + valor);
b = 0.1 * (1 + valor);
l = 0.3 * (1.5 - valor/2);
g = 9.81;

fprintf("\n[Parameters] numero=%d, valor=%.6f\n", numero, valor)
fprintf("  M=%.4f kg  m=%.4f kg  b=%.4f N·s/m  l=%.4f m  g=%.2f m/s^2\n", M, m, b, l, g)

%% Continuous linearized model
A = [0, 1, 0, 0;
     0, -b/M, -m*g/M, 0;
     0, 0, 0, 1;
     0, b/(M*l), (M+m)*g/(M*l), 0];

B_ss = [0; 1/M; 0; -1/(M*l)];

fprintf("\n[Continuous model]\n")
fprintf("  A:\n"); disp(A);
fprintf("  B:\n"); disp(B_ss);

%% Discretization (ZOH)
T = 0.02;
sys_c = ss(A, B_ss, eye(4), zeros(4,1));
sys_d = c2d(sys_c, T, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;

fprintf("[Discrete model] T=%.3f s\n", T)
fprintf("  Ad:\n"); disp(Ad);
fprintf("  Bd:\n"); disp(Bd);

%% Controllability
C_matrix = ctrb(Ad, Bd);
rank_C = rank(C_matrix);
fprintf("[Controllability] rank=%d/4\n", rank_C);
if rank_C < 4, warning("System is NOT controllable!"); end

%% Open-loop eigenvalues
ol_eigs = eig(Ad);
fprintf("[Open-loop] eigenvalues: "); fprintf("%.4f  ", ol_eigs); fprintf("\n");
fprintf("  Unstable modes (|z|>=1): %d\n", sum(abs(ol_eigs) >= 1));

%% DLQR design
Q = diag([100, 80, 10, 1]);
R = 30;

[K, ~, poles_d] = dlqr(Ad, Bd, Q, R);

fprintf("\n[DLQR] Q=diag("); fprintf("%.0f ", diag(Q)); fprintf("), R=%.0f\n", R);
fprintf("  K: "); fprintf("%.4f  ", K); fprintf("\n");
fprintf("  Closed-loop poles: "); fprintf("%.4f  ", poles_d); fprintf("\n");
fprintf("  Max |pole| = %.6f\n", max(abs(poles_d)));

%% Simulation
t_final = 15;
N_steps = round(t_final / T);

x0 = [0; 0; 0; 0];
X_F = [0.5; 0; 0; 0];
U_F = 0;

x_hist = zeros(4, N_steps+1);
u_hist = zeros(1, N_steps);
t_hist = (0:N_steps) * T;

x_hist(:,1) = x0;

ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

for k = 1:N_steps
    u_k = -K * (x_hist(:,k) - X_F) + U_F;
    u_hist(k) = u_k;

    t_start = (k-1) * T;
    t_end = k * T;
    [~, x_ode] = ode45(@(t, s) tarea2_dyn(t, s, u_k, M, m, l, b, g), ...
                        [t_start, t_end], x_hist(:,k), ode_opts);
    x_hist(:, k+1) = x_ode(end, :)';
end

%% Results
x_final = x_hist(:, end);
max_force = max(abs(u_hist));
idx_3s = find(t_hist >= 3.0, 1);
theta_at_3s = x_hist(3, idx_3s);

x1 = x_hist(1, :);
max_x1 = max(x1);
overshoot_val = max_x1 - 0.5;
overshoot_pct = overshoot_val / 0.5 * 100;

settling_band = 0.05 * 0.5;
settled = abs(x1 - 0.5) <= settling_band;
exit_idx = find(~settled);
if isempty(exit_idx)
    settling_time = 0;
else
    settling_time = t_hist(min(exit_idx(end) + 1, length(t_hist)));
end

idx_4s = round(4.0 / T) + 1;
x1_at_4s = x_hist(1, idx_4s);

fprintf("\n[Results - Part b]\n");
fprintf("  Final state:      x1=%.3f  x2=%.3f  x3=%.6f  x4=%.6f\n", ...
    x_final(1), x_final(2), x_final(3), x_final(4));
fprintf("  Max |F|         = %.3f N (limit 1.2)   %s\n", max_force, passfail(max_force < 1.2));
fprintf("  theta at t=3.0s = %.6f rad\n", theta_at_3s);
fprintf("  Max x1          = %.3f\n", max_x1);
fprintf("  Overshoot       = %.4f%% (limit 1%%)    %s\n", overshoot_pct, passfail(overshoot_pct < 1.0));
fprintf("  Settling (5%%)   = %.3f s (limit 4.0)   %s\n", settling_time, passfail(settling_time < 4.0));
fprintf("  x1(k=200) t=4s  = %.3f\n", x1_at_4s);
fprintf("  Weights: r1=%.0f  q1=%.0f  q2=%.0f  q3=%.0f  q4=%.0f\n", R, Q(1,1), Q(2,2), Q(3,3), Q(4,4));

%% Plots
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

figure('Position', [100, 100, 600, 400], 'Color', 'w');
stairs(t_hist(1:end-1), u_hist, 'b-', 'LineWidth', 1.5);
hold on;
yline(1.2, 'r--', 'LineWidth', 1); yline(-1.2, 'r--', 'LineWidth', 1);
hold off;
xlabel('$t$ [s]'); ylabel('$F$ [N]');
title('Control Input $u(k)$');
grid on; box off;

figure('Position', [100, 550, 900, 600], 'Color', 'w');

subplot(2,2,1);
stairs(t_hist, x_hist(1,:), 'b-', 'LineWidth', 1.5);
hold on;
yline(0.5, 'k:', 'LineWidth', 1);
yline(0.505, 'r--', 'LineWidth', 1);
hold off;
xlabel('$t$ [s]'); ylabel('$x$ [m]');
title('Cart Position'); grid on; box off;
legend('$x(t)$', 'Reference', '1\% overshoot', 'Location', 'best');

subplot(2,2,2);
stairs(t_hist, x_hist(2,:), 'b-', 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$\dot{x}$ [m/s]');
title('Cart Velocity'); grid on; box off;

subplot(2,2,3);
stairs(t_hist, x_hist(3,:), 'b-', 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$\theta$ [rad]');
title('Pole Angle'); grid on; box off;

subplot(2,2,4);
stairs(t_hist, x_hist(4,:), 'b-', 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('$\dot{\theta}$ [rad/s]');
title('Pole Angular Velocity'); grid on; box off;

sgtitle('Tarea 2b: Discrete LQR Setpoint Tracking', 'FontSize', 16, 'FontWeight', 'bold');

%% Helper
function s = passfail(cond)
    if cond, s = '[OK]'; else, s = '[VIOLATED]'; end
end

%% Nonlinear dynamics
function dstate = tarea2_dyn(~, state, F, M, m, l, b, g)
    x_dot = state(2);
    theta = state(3);
    theta_dot = state(4);

    ct = cos(theta);
    st = sin(theta);

    M_mat = [(M+m),  m*l*ct;
             m*l*ct, m*l^2];

    rhs = [F - b*x_dot + m*l*st*theta_dot^2;
           m*l*g*st];

    acc = M_mat \ rhs;

    dstate = [x_dot; acc(1); theta_dot; acc(2)];
end
