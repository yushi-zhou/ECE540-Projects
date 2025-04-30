
%% Parameters
lambda = 1;           % Wavelength
k = 2 * pi / lambda;  % Wave number
a_values = [1, 3];    % Side lengths in wavelengths (1λ and 3λ)
gamma = 1.7180724;
e = 2.71828;  % Euler's number

for a_idx = 1:length(a_values)
    a = a_values(a_idx);  % Current side length
    M = ceil(30 * a);     % Segments per side (approx. λ/10 per segment)
    N = 4 * M;            % Total segments
    Delta = a / M;        % Segment length

    x_centers = zeros(N, 1);
    y_centers = zeros(N, 1);
    
    nx = zeros(N, 1);
    ny = zeros(N, 1);
    
    for n = 1:N
        if n <= M
            % Bottom side: x from -a/2 to a/2, y = -a/2
            x_centers(n) = -a/2 + (n - 0.5) * Delta;
            y_centers(n) = -a/2;
            nx(n) = 0;   
            ny(n) = -1;  
        elseif n <= 2*M
            % Right side: x = a/2, y from -a/2 to a/2
            x_centers(n) = a/2;
            y_centers(n) = -a/2 + (n - M - 0.5) * Delta;
            nx(n) = 1;   
            ny(n) = 0;
        elseif n <= 3*M
            % Top side: x from a/2 to -a/2, y = a/2
            x_centers(n) = a/2 - (n - 2*M - 0.5) * Delta;
            y_centers(n) = a/2;
            nx(n) = 0;   
            ny(n) = 1;
        else
            % Left side: x = -a/2, y from a/2 to -a/2
            x_centers(n) = -a/2;
            y_centers(n) = a/2 - (n - 3*M - 0.5) * Delta;
            nx(n) = -1;  
            ny(n) = 0;
        end
    end

    Z = zeros(N, N);
    Z0 = 377;  %impedance of free space (Ohm)
    
    for m = 1:N
        for n = 1:N
            if m == n
                Z(m, n) = -0.5; % Diagonal term for TE is -1/2
            else
                dx = x_centers(m) - x_centers(n);
                dy = y_centers(m) - y_centers(n);
                R = sqrt(dx^2 + dy^2);
                
                if n <= M
                    t_start = -a/2 + (n - 1) * Delta;
                    t_end = -a/2 + n * Delta;
                    t = (t_start + t_end) / 2;
                    Z(m, n) = Delta * gradient_green_function_TE(x_centers(m), y_centers(m), t, -a/2, nx(n), ny(n), k);
                elseif n <= 2*M
                    t_start = -a/2 + (n - M - 1) * Delta;
                    t_end = -a/2 + (n - M) * Delta;
                    t = (t_start + t_end) / 2;
                    Z(m, n) = Delta * gradient_green_function_TE(x_centers(m), y_centers(m), a/2, t, nx(n), ny(n), k);
                elseif n <= 3*M
                    t_start = -a/2 + (n - 2*M - 1) * Delta;
                    t_end = -a/2 + (n - 2*M) * Delta;
                    t = (t_start + t_end) / 2;
                    Z(m, n) = Delta * gradient_green_function_TE(x_centers(m), y_centers(m), (a/2 - (t+a/2)), a/2, nx(n), ny(n), k);
                else
                    t_start = -a/2 + (n - 3*M - 1) * Delta;
                    t_end = -a/2 + (n - 3*M) * Delta;
                    t = (t_start + t_end) / 2;
                    Z(m, n) = Delta * gradient_green_function_TE(x_centers(m), y_centers(m), -a/2, (a/2 - (t+a/2)), nx(n), ny(n), k);
                end
                
                %Z(m, n) = integral(integrand, t_start, t_end, 'RelTol', 1e-3);
            end
        end
    end

    phi_i = 0;  % Incidence angle
    V = zeros(N, 1);  % Incident field vector
    for n = 1:N
        H_inc = exp(-1j * k * (x_centers(n) * cos(phi_i) + y_centers(n) * sin(phi_i)));
                
        V(n) = H_inc;
    end

    J = Z \ V; 

    %%calculate total field at anywhere
    grid_size = 100;
    x_grid = linspace(-5*a, 5*a, grid_size);
    y_grid = linspace(-5*a, 5*a, grid_size);
    [X, Y] = meshgrid(x_grid, y_grid);
    Hz_total = zeros(size(X));
    Hz_st = zeros(size(X));

    for ii = 1:numel(X)
        x = X(ii);
        y = Y(ii);
        if (x >= -a/2 && x <= a/2 && y >= -a/2 && y <= a/2)
            Hz_total(ii) = 0;  
        else
            % Scattered field
            Hz_s = 0;
            for n = 1:N
                R = sqrt((x - x_centers(n))^2 + (y - y_centers(n))^2);
                if R < 1e-10
                    R = 1e-10;  
                end

                dx = x - x_centers(n);
                dy = y - y_centers(n);

                dG_dn = k/(4*1j) * besselh(1, 2, k*R) * (dx*nx(n) + dy*ny(n))/R;

                Hz_s = Hz_s - J(n) * dG_dn * Delta;
            end
            Hz_i = exp(-1j * k * (x * cos(phi_i) + y * sin(phi_i)));
            Hz_total(ii) = Hz_i + Hz_s;
            Hz_st(ii) = Hz_s;
        end
    end

    %% Plot total field magnitude
    figure;
    mask = (X >= -a/2 & X <= a/2 & Y >= -a/2 & Y <= a/2);
    
    Hz_plot = abs(Hz_total);
    
    pcolor(X, Y, Hz_plot);
    shading interp;
    
    colormap_with_white = colormap;
    hold on;
    
    rectangle('Position', [-a/2, -a/2, a, a], 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Apply white mask to interior points
    white_patch = patch([-a/2 a/2 a/2 -a/2], [-a/2 -a/2 a/2 a/2], 'w');
    set(white_patch, 'EdgeColor', 'none');
    
    colorbar;
    title(['Total Field Magnitude (TE), a = ', num2str(a), '\lambda']);
    xlabel('x (\lambda)');
    ylabel('y (\lambda)');
    axis equal;

    %% bistatic scattering width
    theta = linspace(0, 360, 361);  % Scattering angles in degrees
    f_theta = zeros(size(theta));
    for tt = 1:length(theta)
        th = theta(tt) * pi / 180;  % Convert to radians
        sum_term = 0;
        for n = 1:N
            sum_term = sum_term + J(n) * exp(1j * k * (x_centers(n) * cos(th) + y_centers(n) * sin(th)));
        end
        f_theta(tt) = Delta * sum_term;
    end

    sigma_2D = 2*pi * abs(f_theta).^2; % Proper 2D scattering width formula
    sigma_2D_db = 10 * log10(sigma_2D / lambda); % Normalized by wavelength

    figure;
    plot(theta, sigma_2D_db);
    xlabel('Scattering Angle (degrees)');
    xlim([0, 360]);
    ylabel('Bistatic Scattering Width \sigma_{2D}/\lambda (dB)');
    title(['Bistatic Scattering Width vs. Angle (TE), a = ', num2str(a), '\lambda']);
    grid on;

    %% Surface Current Density
    figure;
    perimeter_pos = zeros(N, 1);
    for n = 1:N
        if n <= M
            perimeter_pos(n) = (n - 1) * Delta;
        elseif n <= 2*M
            perimeter_pos(n) = M*Delta + (n - M - 1) * Delta;
        elseif n <= 3*M
            perimeter_pos(n) = 2*M*Delta + (n - 2*M - 1) * Delta;
        else
            perimeter_pos(n) = 3*M*Delta + (n - 3*M - 1) * Delta;
        end
    end
    
    Js_t = abs(J); 

    plot(perimeter_pos, Js_t, 'LineWidth', 1.5);
    xlabel('Position along boundary (\lambda)');
    ylabel('|J_{s,t}|');
    title(['Induced Surface Current Density (TE), a = ', num2str(a), '\lambda']);
    grid on;
    xticks([0 M*Delta 2*M*Delta 3*M*Delta 4*M*Delta]);
    xticklabels({'0', 'a', '2a', '3a', '4a'});
end


%Green's function gradient calculation
function val = gradient_green_function_TE(x_obs, y_obs, x_src, y_src, nx, ny, k)
    dx = x_obs - x_src;
    dy = y_obs - y_src;
    R = sqrt(dx.^2 + dy.^2);  % Use element-wise power .^ instead of ^
    
    if any(R < 1e-10)
        R(R < 1e-10) = 1e-10;  % Replace very small values with a minimum
    end
    
    dH_dr = k/(4*1j) * besselh(1, 2, k*R);
    
    val = dH_dr .* (dx*nx + dy*ny) ./ R;
end