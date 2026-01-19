function x_next = rk4_step(f, t, x, dt)
    % RK4_STEP Performs one step of 4th-order Runge-Kutta integration
    %
    % Inputs:
    %   f  - Function handle for derivative: dx/dt = f(t, x)
    %   t  - Current time
    %   x  - Current state vector
    %   dt - Time step
    %
    % Output:
    %   x_next - State at next time step

    k1 = f(t, x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt, x + dt * k3);

    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end
