function [t, y, x_hat, sigma] = cart_pole_l1_solver(f, params, x0, sim_time, r)
    % CART_POLE_L1_SOLVER L1 adaptive control solver for cart-pole systems
    %
    % Inputs:
    %   f        - Dynamics function handle
    %   params   - Parameter structure
    %   x0       - Initial state
    %   sim_time - Simulation time
    %   r        - Reference function handle r(t) (default: @(t) 0)

    arguments
        f
        params
        x0
        sim_time
        r = @(t) 0
    end

    % Initialize
    x = x0;                       % Plant state
    x_hat = x0;                   % Predictor state
    u_filt = 0;                   % Filtered control

    % Determine number of states and unmatched uncertainties
    n_states = length(x0);
    n_um = size(params.Bum, 2);   % Number of unmatched uncertainties

    sigma_hat_m = 0;              % Matched estimate (scalar)
    sigma_hat_um = zeros(n_um, 1); % Unmatched estimates (vector)

    % Time vector
    t = 0:params.Ts:sim_time;
    n_steps = length(t);

    % Preallocate
    y = zeros(n_steps, n_states);
    x_hat_hist = zeros(n_steps, n_states);
    sigma_hist = zeros(n_steps, 1 + n_um);  % [sigma_m, sigma_um(:)']
    
    % Check for mass step change parameters
    has_single_mass_step = isfield(params, 'm_step_time') && isfield(params, 'm_step_value');
    has_double_mass_step = isfield(params, 'm1_step_time') && isfield(params, 'm1_step_value');
    mass_changed = false;

    % Simulation loop
    for i = 1:n_steps
        % Apply mass step change if needed
        if ~mass_changed
            if has_single_mass_step && t(i) >= params.m_step_time
                params.m = params.m_step_value;
                mass_changed = true;
                fprintf('Mass step applied at t=%.2f: m = %.3f kg\n', t(i), params.m);
            elseif has_double_mass_step && t(i) >= params.m1_step_time
                params.m1 = params.m1_step_value;
                mass_changed = true;
                fprintf('Mass step applied at t=%.2f: m1 = %.3f kg\n', t(i), params.m1);
            end
        end

        % Store
        y(i,:) = x';
        x_hat_hist(i,:) = x_hat';
        sigma_hist(i,:) = [sigma_hat_m; sigma_hat_um]';

        % State error
        x_tilde = x_hat - x;

        % Compute sigma matched and sigma unmatched estimation
        sigma_update = params.gamma * x_tilde;
        sigma_hat_m = sigma_update(1);
        sigma_hat_um = sigma_update(2:end);

        % Adaptive control law (feedforward + uncertainty compensation)
        u_ad_raw = params.k_g * r(t(i)) - sigma_hat_m;

        % Discrete low-pass filter
        a = exp(params.Ts * params.wf);
        u_filt = (1 - a) * u_filt + a * u_ad_raw;

        % Total control to plant (baseline feedback + filtered adaptive)
        u_total = -params.K_lqr * x_hat + u_filt;

        % Solve plant dynamics using RK4
        f_plant = @(t_val, x_val) f(t_val, x_val, params, @(t,x) u_total);
        x = rk4_step(f_plant, t(i), x, params.Ts);

        % Predict the state using RK4
        f_pred = @(t_val, x_hat_in) params.Am * x_hat_in + params.Bm * (u_filt + sigma_hat_m) + params.Bum * sigma_hat_um;
        x_hat = rk4_step(f_pred, t(i), x_hat, params.Ts);
    end
    
    x_hat = x_hat_hist;
    sigma = sigma_hist;
end
