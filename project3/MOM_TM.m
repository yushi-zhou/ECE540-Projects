% MATLAB code to solve scattering from a square conducting cylinder using MOM
% TM polarization, 2D cross-section problem
% Plots field pattern in the domain and amplitude vs. scattering angle

% Parameters
lambda = 1;           % Wavelength
k = 2 * pi / lambda;  % Wave number
a_values = [1, 3];    % Side lengths in wavelengths (1λ and 3λ)
gamma = 1.7180724;
e = 2.71828;  % Euler's number

% Discretization levels for convergence study
M_factors = [10, 20, 30, 50, 100]; 

for a_idx = 1:length(a_values)
    a = a_values(a_idx);  % Current side length
    
    % For convergence study
    Js_all = cell(length(M_factors), 1);  % Store current densities for all discretizations
    perimeter_pos_all = cell(length(M_factors), 1);  % Store perimeter positions
    
    % For main analysis with M = ceil(30*a)
    M_main_idx = find(M_factors == 30);
    
    % Loop over different discretization levels
    for m_idx = 1:length(M_factors)
        M = ceil(M_factors(m_idx) * a);  % Number of segments per side
        N = 4 * M;                     % Total segments
        Delta = a / M;                 % Segment length
        
        % Define segment centers along the square perimeter (counterclockwise)
        x_centers = zeros(N, 1);
        y_centers = zeros(N, 1);
        
        
        for n = 1:N
            if n <= M
                % Bottom side: x from -a/2 to a/2, y = -a/2
                x_centers(n) = -a/2 + (n - 0.5) * Delta;
                y_centers(n) = -a/2;
            elseif n <= 2*M
                % Right side: x = a/2, y from -a/2 to a/2
                x_centers(n) = a/2;
                y_centers(n) = -a/2 + (n - M - 0.5) * Delta;
            elseif n <= 3*M
                % Top side: x from a/2 to -a/2, y = a/2
                x_centers(n) = a/2 - (n - 2*M - 0.5) * Delta;
                y_centers(n) = a/2;
            else
                % Left side: x = -a/2, y from a/2 to -a/2
                x_centers(n) = -a/2;
                y_centers(n) = a/2 - (n - 3*M - 0.5) * Delta;
            end
        end

        Z = zeros(N, N);
        Z0 = 377;  % impedance of free space (Ohm)
        factor = k*Z0/4;
        for m = 1:N
            for n = 1:N
                if m == n
                    Z(m, n) = factor * Delta * (1 - 1j * 2 / pi * log(k * Delta * gamma / (4*e)));
                else
                    if n <= M
                        t_start = -a/2 + (n - 1) * Delta;
                        t_end = -a/2 + n * Delta;
                        t = (t_start + t_end) / 2;
                        Z(m, n) = Delta * factor * besselh(0, 2, k * sqrt((x_centers(m) - t).^2 + (y_centers(m) - (-a/2)).^2));
                    elseif n <= 2*M
                        t_start = -a/2 + (n - M - 1) * Delta;
                        t_end = -a/2 + (n - M) * Delta;
                        t = (t_start + t_end) / 2;
                        Z(m, n) = Delta * factor * besselh(0, 2, k * sqrt((x_centers(m) - (a/2)).^2 + (y_centers(m) - t).^2));
                    elseif n <= 3*M
                        t_start = -a/2 + (n - 2*M - 1) * Delta;
                        t_end = -a/2 + (n - 2*M) * Delta;
                        t = (t_start + t_end) / 2;
                        Z(m, n) = Delta * factor * besselh(0, 2, k * sqrt((x_centers(m) - (a/2 - (t+a/2))).^2 + (y_centers(m) - (a/2)).^2));
                    else
                        t_start = -a/2 + (n - 3*M - 1) * Delta;
                        t_end = -a/2 + (n - 3*M) * Delta;
                        t = (t_start + t_end) / 2;
                        Z(m, n) = Delta * factor * besselh(0, 2, k * sqrt((x_centers(m) - (-a/2)).^2 + (y_centers(m) - (a/2 - (t+a/2))).^2));
                    end
                end
            end
        end

        phi_i = 0;  % Incidence angle 
        V = zeros(N, 1);  
        for n = 1:N
            V(n) = exp(-1j * k * (x_centers(n) * cos(phi_i) + y_centers(n) * sin(phi_i)));
        end

        I = Z \ V;  

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
        
        Js_z = abs(I);  
        
        perimeter_pos_normalized = perimeter_pos / (4*Delta*M);
        
        Js_all{m_idx} = Js_z;
        perimeter_pos_all{m_idx} = perimeter_pos_normalized * 4*a;  
        
        % Only do the full analysis for M = ceil(30*a)
        if M_factors(m_idx) == 30
            %%calculate total field at anywhere
            grid_size = 100;
            x_grid = linspace(-5*a, 5*a, grid_size);
            y_grid = linspace(-5*a, 5*a, grid_size);
            [X, Y] = meshgrid(x_grid, y_grid);
            E_total = zeros(size(X));
            E_st = zeros(size(X));

            for ii = 1:numel(X)
                x = X(ii);
                y = Y(ii);
                if (x >= -a/2 && x <= a/2 && y >= -a/2 && y <= a/2)
                    E_total(ii) = 0;  % Field inside conductor is zero
                else
                    % Scattered field
                    E_s = 0;
                    for n = 1:N
                        R = sqrt((x - x_centers(n))^2 + (y - y_centers(n))^2);
                        E_s = E_s + I(n) * (-k*Z0/4 * besselh(0, 2, k * R)) * Delta;
                    end
                    % Incident field
                    E_i = exp(-1j * k * (x * cos(phi_i) + y * sin(phi_i)));

                    E_total(ii) = E_i + E_s;
                    E_st(ii) = E_s;
                end
            end

            %% --- Plot Geometry and Surface Discretization ---
            figure;
            hold on;
            rectangle('Position', [-a/2, -a/2, a, a], 'EdgeColor', 'k', 'LineWidth', 1.5);
            % Plot the segment centers
            plot(x_centers, y_centers, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
            axis equal;
            xlabel('x (\lambda)');
            ylabel('y (\lambda)');
            title(['Geometry and Surface Discretization, a = ', num2str(a), '\lambda']);
            legend('Segment Centers', 'Normal Vectors');
            grid on;

            %% Plot total field magnitude
            figure;
            mask = (X >= -a/2 & X <= a/2 & Y >= -a/2 & Y <= a/2);
            
            E_plot = abs(E_total);
            
            pcolor(X, Y, E_plot);
            shading interp;
            
            hold on;
            rectangle('Position', [-a/2, -a/2, a, a], 'EdgeColor', 'k', 'LineWidth', 1.5);
            
            % Apply white mask to interior points
            white_patch = patch([-a/2 a/2 a/2 -a/2], [-a/2 -a/2 a/2 a/2], 'w');
            set(white_patch, 'EdgeColor', 'none');
            
            colorbar;
            title(['Total Field Magnitude (TM), a = ', num2str(a), '\lambda']);
            xlabel('x (\lambda)');
            ylabel('y (\lambda)');
            axis equal;

            %% Bistatic Scattering Width ---
            theta = linspace(0, 360, 361);  % Scattering angles in degrees
            sigma_2D = zeros(size(theta));  % Bistatic scattering width

            for tt = 1:length(theta)
                th = theta(tt) * pi / 180;  % Convert to radians
                sum_term = 0;
                for n = 1:N
                    sum_term = sum_term + I(n) * exp(1j * k * (x_centers(n) * cos(th) + y_centers(n) * sin(th)));
                end
                f_theta = Delta * sum_term;
                sigma_2D(tt) = (abs(f_theta)^2) / (2 * pi);  % Bistatic scattering width formula
            end

            sigma_2D_db = 10 * log10(sigma_2D);

            figure;
            plot(theta, sigma_2D_db, 'LineWidth', 1.5);
            xlabel('Scattering Angle (degrees)');
            xlim([0, 360]);
            ylabel('\sigma_{2D}/\lambda (dB)');
            title(['Bistatic Scattering Width (TM), a = ', num2str(a), '\lambda']);
            grid on;

            %% --- Plot Induced Surface Current Density (single M value) ---
            figure;
            
            plot(perimeter_pos, Js_z, 'LineWidth', 1.5);
            xlabel('Position along boundary (\lambda)');
            ylabel('|J_{s,z}|');
            title(['Induced Surface Current Density (TM), a = ', num2str(a), '\lambda, M = ', num2str(M_factors(m_idx)), '*a']);
            grid on;
            xticks([0 M*Delta 2*M*Delta 3*M*Delta 4*M*Delta]);
            xticklabels({'0', 'a', '2a', '3a', '4a'});
        end
    end
    
    %% --- Plot Convergence of Current Density for Different M Values ---
    figure;
    hold on;
    
    % Define colors and line styles for the convergence plot
    colors = {'b', 'r', 'g', 'm', 'black'};
    line_styles = {'-', '-', '-', '-', '-'};
    
    legend_entries = cell(length(M_factors), 1);
    
    for m_idx = 1:length(M_factors)
        plot(perimeter_pos_all{m_idx}, Js_all{m_idx}, [colors{m_idx}, line_styles{m_idx}], 'LineWidth', 1.5);
        legend_entries{m_idx} = ['M = ', num2str(M_factors(m_idx)), '*a'];
    end
    
    xlabel('Position along boundary (\lambda)');
    ylabel('|J_{s,z}|');
    title(['Convergence of Induced Surface Current Density (TM), a = ', num2str(a), '\lambda']);
    legend(legend_entries, 'Location', 'best');
    grid on;
    
    % Add x-ticks at corner positions
    xticks([0 a 2*a 3*a 4*a]);
    xticklabels({'0', 'a', '2a', '3a', '4a'});
end