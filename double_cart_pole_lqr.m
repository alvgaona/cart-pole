clear; clc; close all;

%% Parameters of the model
M = 1;      % Mass of the cart [kg]
m1 = 1;     % Mass of the first pole [kg]
m2 = 1;     % Mass of the second pole [kg]
l1 = 1;     % Length of the first pole [m]
l2 = 1;     % Length of the second pole [m]
b = 0.01;   % Friction coefficiente of the cart
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

