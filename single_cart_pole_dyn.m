function dstate = single_cart_pole_dyn(t, state, params, control_func, dist_func)
% SINGLE_CART_POLE DYN Computes the dynamics of single pendulum cart pole
%
% Inputs:
%   t - current time
%   state - state vector [x; theta; dx; dtheta]
%           x: forward position, theta: pole angle, dtheta: pole angular velocity
%   params - structure containing robot parameters
%   control_func - function handle for control input: u = control_func(t, state)
%   dist_func - function handle for disturbances: [F_cart; tau_pole] = dist_func(t, state)
%               F_cart: force on cart [N], tau_pole: torque on pole [N·m]
%               (default: @(t,x) [0; 0])
%
% Output:
%   dstate - time derivative of state vector

arguments
    t
    state
    params
    control_func
    dist_func = @(t, x) [0; 0]
end

% Extract state variables
x = state(1);
theta = state(2);
dx = state(3);
dtheta = state(4);

% Extract parameters
M = params.M;
m = params.m;
l = params.l;
b = params.b;
g = params.g;

% Mass matrix
a11 = M+m;
a12 = -m*l*cos(theta);
a21 = a12;
a22 = m*l^2;

M = [
    a11 a12;
    a21 a22;
];

% Coriolis matrix
c12 = m*l*sin(theta)*dtheta;

C = [
    0 c12;
    0 0
];

% Damping matrix
D = [
    b 0;
    0 0
];

% Input matrix
B = [1; 0];

% Gravity vector
G = [0; -m*l*g*sin(theta)];

% Get control input
u = control_func(t, state);

% Get disturbances: [F_cart; tau_pole]
d = dist_func(t, state);
F_cart = d(1);    % Force on cart [N]
tau_pole = d(2);  % Torque on pole [N·m]

% Disturbance input matrix
B_dist = [1; 0; 1/l];  % F_cart affects cart, tau_pole affects pole rotation

% Velocity vector
qdot = [dx; dtheta];

% Equation of motion: M*qddot + C*qdot + D*qdot + G = B*u + B_dist*[F_cart; tau_pole]
% Solve for qddot
qddot = M \ (B*u + [F_cart; tau_pole] - C*qdot - D*qdot - G);

% Construct derivative of state vector
dstate = [
    dx;           % d/dt(x) = dx
    dtheta;       % d/dt(theta) = dtheta
    qddot(1);     % d/dt(dx) = ddx
    qddot(2);     % d/dt(dtheta) = ddtheta
];

end

