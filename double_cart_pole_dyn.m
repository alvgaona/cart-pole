function dstate = double_cart_pole_dyn(t, state, params, control_func, dist_func)
% DOUBLE_CART_POLE DYN Computes the dynamics of double pendulum cart pole
%
% Inputs:
%   t - current time
%   state - state vector [x; theta1; theta2; dx; dtheta1; dtheta2]
%   params - structure containing robot parameters
%   control_func - function handle for control input: u = control_func(t, state)
%   dist_func - function handle for disturbances: [F_cart; tau_pole1; tau_pole2] = dist_func(t, state)
%               F_cart: force on cart [N], tau_pole1/2: torques on poles [N·m]
%               (default: @(t,x) [0; 0; 0])
%
% Output:
%   dstate - time derivative of state vector

arguments
    t
    state
    params
    control_func
    dist_func = @(t, x) [0; 0; 0]
end

% Extract state variables
x = state(1);
theta1 = state(2);
theta2 = state(3);
dx = state(4);
dtheta1 = state(5);
dtheta2 = state(6);

% Extract parameters
M = params.M;
m1 = params.m1;
m2 = params.m2;
l1 = params.l1;
l2 = params.l2;
b = params.b;
g = params.g;

% Mass matrix
a11 = M+m1+m2;
a12 = -(m1+m2)*l1*cos(theta1);
a13 = -m2*l2*cos(theta2);
a21 = a12;
a22 = (m1+m2)*l1^2;
a23 = m2*l1*l2*cos(theta1-theta2);
a31 = a13;
a32 = a23;
a33 = m2*l2^2;

M = [
    a11 a12 a13;
    a21 a22 a23;
    a31 a32 a33
];

cond_num = cond(M);
if cond_num > 1e10
    warning('Mass matrix ill-conditioned (cond=%.2e) at t=%.3f', cond_num, t);
end

% Coriolis matrix
c12 = (m1+m2)*l1*sin(theta1)*dtheta1;
c13 = m2*l2*sin(theta2)*dtheta2;
c23 = m2*l1*l2*sin(theta1-theta2)*dtheta2;
c32 = -m2*l1*l2*sin(theta1-theta2)*dtheta1;

C = [
    0 c12 c13;
    0 0 c23;
    0 c32 0
];

% Damping matrix
D = [
    b 0 0;
    zeros(2,3)
];

% Input matrix
B = [1; 0; 0];

% Gravity vector
G = [0; -(m1+m2)*l1*g*sin(theta1); -m2*l2*g*sin(theta2)];

% Get control input
u = control_func(t, state);

% Get disturbances: [F_cart; tau_pole1; tau_pole2]
d = dist_func(t, state);
F_cart = d(1);      % Force on cart [N]
tau_pole1 = d(2);   % Torque on pole 1 [N·m]
tau_pole2 = d(3);   % Torque on pole 2 [N·m]

% Velocity vector
qdot = [dx; dtheta1; dtheta2];

% Equation of motion: M*qddot + C*qdot + D*qdot + G = B*u + [F_cart; tau_pole1; tau_pole2]
% Solve for qddot
qddot = M \ (B*u + [F_cart; tau_pole1; tau_pole2] - C*qdot - D*qdot - G);

% Construct derivative of state vector
dstate = [
    dx;           % d/dt(x) = dx
    dtheta1;      % d/dt(theta1) = dtheta1
    dtheta2;      % d/dt(theta2) = dtheta2
    qddot(1);     % d/dt(dx) = ddx
    qddot(2);     % d/dt(dtheta1) = ddtheta1
    qddot(3);     % d/dt(dtheta2) = ddtheta2
];

end

