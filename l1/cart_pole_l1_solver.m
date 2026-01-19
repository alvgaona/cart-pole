function [t, y, x_hat, sigma] = cart_pole_l1_solver(f, params, x0, sim_time, r, dist_func)
    % CART_POLE_L1_SOLVER L1 adaptive control solver for cart-pole systems
    %
    % Inputs:
    %   f         - Dynamics function handle
    %   params    - Parameter structure
    %   x0        - Initial state
    %   sim_time  - Simulation time
    %   r         - Reference function handle r(t) (default: @(t) 0)
    %   dist_func - Disturbance function handle (default: @(t,x) zeros)

    arguments
        f
        params
        x0
        sim_time
        r = @(t) 0
        dist_func = []
    end

    % Set default disturbance function if not provided
    if isempty(dist_func)
        n_states = length(x0);
        if n_states == 4
            dist_func = @(t, x) [0; 0];  % Single cart-pole: [F_cart; tau_pole]
        else
            dist_func = @(t, x) [0; 0; 0];  % Double cart-pole: [F_cart; tau_pole1; tau_pole2]
        end
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
    
    % Check for parameter step changes
    has_single_mass_step = isfield(params, 'm_step_time') && isfield(params, 'm_step_value');
    has_single_length_step = isfield(params, 'l_step_time') && isfield(params, 'l_step_value');
    has_m1_step = isfield(params, 'm1_step_time') && isfield(params, 'm1_step_value');
    has_m2_step = isfield(params, 'm2_step_time') && isfield(params, 'm2_step_value');
    has_l1_step = isfield(params, 'l1_step_time') && isfield(params, 'l1_step_value');
    has_l2_step = isfield(params, 'l2_step_time') && isfield(params, 'l2_step_value');
    has_friction_step = isfield(params, 'b_step_time') && isfield(params, 'b_step_value');

    mass_changed = false;
    length_changed = false;
    m1_changed = false;
    m2_changed = false;
    l1_changed = false;
    l2_changed = false;
    friction_changed = false;

    % Simulation loop
    for i = 1:n_steps
        % Apply single pole mass step change if needed
        if ~mass_changed && has_single_mass_step && t(i) >= params.m_step_time
            params.m = params.m_step_value;
            mass_changed = true;
            fprintf('Mass step applied at t=%.2f: m = %.3f kg\n', t(i), params.m);
        end

        % Apply single pole length step change if needed
        if ~length_changed && has_single_length_step && t(i) >= params.l_step_time
            params.l = params.l_step_value;
            length_changed = true;
            fprintf('Length step applied at t=%.2f: l = %.3f m\n', t(i), params.l);
        end

        % Apply m1 step change if needed
        if ~m1_changed && has_m1_step && t(i) >= params.m1_step_time
            params.m1 = params.m1_step_value;
            m1_changed = true;
            fprintf('Mass step applied at t=%.2f: m1 = %.3f kg\n', t(i), params.m1);
        end

        % Apply m2 step change if needed
        if ~m2_changed && has_m2_step && t(i) >= params.m2_step_time
            params.m2 = params.m2_step_value;
            m2_changed = true;
            fprintf('Mass step applied at t=%.2f: m2 = %.3f kg\n', t(i), params.m2);
        end

        % Apply l1 step change if needed
        if ~l1_changed && has_l1_step && t(i) >= params.l1_step_time
            params.l1 = params.l1_step_value;
            l1_changed = true;
            fprintf('Length step applied at t=%.2f: l1 = %.3f m\n', t(i), params.l1);
        end

        % Apply l2 step change if needed
        if ~l2_changed && has_l2_step && t(i) >= params.l2_step_time
            params.l2 = params.l2_step_value;
            l2_changed = true;
            fprintf('Length step applied at t=%.2f: l2 = %.3f m\n', t(i), params.l2);
        end

        % Apply friction step change if needed
        if ~friction_changed && has_friction_step && t(i) >= params.b_step_time
            params.b = params.b_step_value;
            friction_changed = true;
            fprintf('Friction step applied at t=%.2f: b = %.3f N·s/m\n', t(i), params.b);
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
        u_ad = params.k_g * r(t(i)) - sigma_hat_m;

        % Discrete low-pass filter
        a = exp(-params.Ts * params.wf);
        u_filt = a * u_filt + (1 - a) * u_ad;

        % Total control to plant (baseline feedback + filtered adaptive)
        u_total = -params.K_lqr * x + u_filt;

        % Solve plant dynamics using RK4
        f_plant = @(t_val, x_val) f(t_val, x_val, params, @(t,x) u_total, dist_func);
        x = rk4_step(f_plant, t(i), x, params.Ts);

        % Predict the state using RK4
        f_pred = @(t_val, x_hat_in) params.Am * x_hat_in + params.Bm * (u_filt + sigma_hat_m) + params.Bum * sigma_hat_um;
        x_hat = rk4_step(f_pred, t(i), x_hat, params.Ts);
    end
    
    x_hat = x_hat_hist;
    sigma = sigma_hist;
end
