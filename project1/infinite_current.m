% A Matlab program using FDTD to calculate the radiation of a infinite current line, with Yee's grid.
% Author: Yushi Zhou  
function infinite_current()
    %define physical constants
    c = 3e8; %speed of light
    mu0 = 4*pi*1e-7; %permeability of free space
    epsilon0 = 8.85e-12; %permittivity of free space
    freq = 1e9;
    omega = 2*pi*freq;
    lambda = c/freq;
    x_source = 1;
    y_source = 1;

    % Define the grid
    dx = 0.01;
    dy = dx;
    dt = dx/10e8; %  dt<=\frac{dx}{c\sqrt{2}}

    %define range
    x = 0:dx:2;
    y = 0:dy:2;
    %define time range
    t = 0:dt:10/freq;

    %define permitivity and conductivity
    epsilon = ones(length(x),length(y))*epsilon0;
    sigma = zeros(length(x),length(y));

    %define a sin current source, period = 100*dt, (1.5,1)
    J_z = zeros(length(t),length(x),length(y));

    for n=1:length(t)
        J_z(n,round(x_source/dx),round(y_source/dy)) = sin(omega*n*dt)/(dx*dy);
    end

    %define perfect conductor sheet at x=x_conductor
    idx_conductor_x = 50;
    idx_slit_y = 0.3*length(y);

    %define PML 
    PML = round((lambda)/dx);
    sigma_max = -log(0.001)*4/(2*377*PML*dx);
    for i = 1:PML
        for j = 1:length(y)
            sigma(i, j) = sigma_max * ((PML - i) / PML)^3;
            sigma(length(x) - i + 1, j) = sigma_max * ((PML - i) / PML)^3;
        end
    end
    for i = 1:length(x)
        for j = 1:PML
            sigma(i, j) = sigma_max * ((PML - j) / PML)^3;
            sigma(i, length(y) - j + 1) = sigma_max * ((PML - j) / PML)^3;
        end
    end

    % Define PML at the corners
    for i = 1:PML
        for j = 1:PML
            sigma_x = sigma_max * ((PML - i + 0.5) / PML)^5;
            sigma_y = sigma_max * ((PML - j + 0.5) / PML)^5;
            
            sigma(i, j) = max(sigma_x, sigma_y);
            sigma(i, length(y) - j + 1) = max(sigma_x, sigma_y);
            sigma(length(x) - i + 1, j) = max(sigma_x, sigma_y);
            sigma(length(x) - i + 1, length(y) - j + 1) = max(sigma_x, sigma_y);
        end
    end

    %define coefficient alpha, beta
    alpha = epsilon./dt - sigma/2;
    beta = epsilon./dt + sigma/2;

    %define E_z, H_x, H_y on Yee's grid
    E_z = zeros(length(t),length(x),length(y));
    H_x = zeros(length(t),length(x),length(y));
    H_y = zeros(length(t),length(x),length(y));

    %FDTD 
    for n=1:length(t)-1
        %H_x, H_y, from t=3/2 to t=N+1/2
        for i=1:length(x)-1
            for j=1:length(y)-1
                if i>PML && i<length(x)-PML && j>PML && j<length(y)-PML
                    H_x(n+1,i,j+1) = H_x(n,i,j+1) - (E_z(n,i,j+1)-E_z(n,i,j))*dt/(mu0*dy);
                    H_y(n+1,i+1,j) = H_y(n,i+1,j) + (E_z(n,i+1,j)-E_z(n,i,j))*dt/(mu0*dx);
                else %Jin's book 2nd:8.5.62-8.5.65
                    H_x(n+1,i,j+1) = 1./beta(i,j+1)*(alpha(i,j+1).*H_x(n,i,j+1) - epsilon(i,j+1)/mu0/dy*(E_z(n,i,j+1)-E_z(n,i,j)));
                    H_y(n+1,i+1,j) = 1./beta(i+1,j)*(alpha(i+1,j).*H_y(n,i+1,j) + epsilon(i+1,j)/mu0/dx*(E_z(n,i+1,j)-E_z(n,i,j)));
                end
            end
        end

        %E_z, from t=2 to t=N(Known E_z[1])
        for i=1:length(x)-1
            for j=1:length(y)-1
                if i>PML && i<length(x)-PML && j>PML && j<length(y)-PML
                    E_z(n+1,i,j) = 1./beta(i,j)*(alpha(i,j).*E_z(n,i,j) + 1/dx*(H_y(n+1,i+1,j)-H_y(n+1,i,j))-1/dy*(H_x(n+1,i,j+1)-H_x(n+1,i,j))-J_z(n+1,i,j));

                    %add a conductive sheet
                    % if i==idx_conductor_x && (j < idx_slit_y-5 || (j > idx_slit_y+5) && (j < length(y)-idx_slit_y-5) || j > length(y)-idx_slit_y+5)
                    %     E_z(n+1,i,j) = 0;
                    % end
                else
                    E_z(n+1,i,j) = 1./beta(i,j)*(alpha(i,j).*E_z(n,i,j) + 1/dx*(H_y(n+1,i+1,j)-H_y(n+1,i,j))-1/dy*(H_x(n+1,i,j+1)-H_x(n+1,i,j)));

                    %add a conductive sheet
                    % if i==idx_conductor_x && (j < idx_slit_y-5 || (j > idx_slit_y+5) && (j < length(y)-idx_slit_y-5) || j > length(y)-idx_slit_y+5)
                    %     E_z(n+1,i,j) = 0;
                    % end
                end
            end
             
        end
        

    end

    % analytical solution
    E_z_analytical = zeros(size(E_z));
    for n = 1:length(t)
        for i = PML+1:length(x)-PML-1
            for j = PML+1:length(y)-PML-1
                r = sqrt((x(i)-x(x_source/dx))^2 + (y(j)-y(y_source/dy))^2);
                E_z_analytical(n, i, j) = analytical_solution(r, n*dt, omega, mu0, epsilon0, dx, dy);
            end
        end
    end
    error = abs(E_z - E_z_analytical); %error at the beginning is large because system is not in steady yet
    
    % Add this after calculating 'error'
    error_metrics = calculate_error_metrics(E_z, E_z_analytical, PML);
    
    % Analyze PML reflections
    analyze_pml_reflections(E_z, x, y, t, PML);
    
    % plot
    visualize_field(E_z, x, y, t, PML);
 
    save('infinite_current.mat','E_z','H_x','H_y','J_z','x','y','t');

end

function E_z_analytical = analytical_solution(r, t, omega, mu0, epsilon0, dx, dy)
    %Helmholtz solution
    k = omega * sqrt(mu0 * epsilon0);
    E_z_analytical = (mu0*omega/4)*besselh(0, 1, -k*r )*exp(1i*omega*t); % Hankel function
    E_z_analytical = imag(E_z_analytical); 
end
    


function visualize_field(E_z, x, y, t, PML)
    figure;
    hFig = gcf;
    hAx = axes('Parent', hFig);
    hImg = imagesc(x, y, squeeze(E_z(1, :, :))', 'Parent', hAx);
    colorbar;
    colormap(winter);
    %clim([-0.2, 0.2]); % Set color scale limits for better contrast
    title('E_z at different times');
    xlabel('x');
    ylabel('y');
    axis equal;
    axis tight;
    
    % Mark PML boundaries with red lines
    hold on;
    % Draw rectangle to show non-PML region
    non_pml_x = [x(PML+1), x(end-PML), x(end-PML), x(PML+1), x(PML+1)];
    non_pml_y = [y(PML+1), y(PML+1), y(end-PML), y(end-PML), y(PML+1)];
    plot(non_pml_x, non_pml_y, 'r--', 'LineWidth', 2);
    text(x(PML+5), y(PML+5), 'PML Boundary', 'Color', 'r', 'FontWeight', 'bold');
    hold off;

    % Create a slider bar
    hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(t), 'Value', 1, ...
                        'Units', 'normalized', 'Position', [0.2 0.01 0.6 0.05], ...
                        'Callback', @(src, event) updateImage());

    % update the image with the bar
    addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateImage());
    function updateImage()
        n = round(get(hSlider, 'Value'));
        set(hImg, 'CData', squeeze(E_z(n, :, :))');
        title(hAx, sprintf('E_z at t = %.2f', t(n)));
        
        % Ensure PML boundary remains visible
        hold on;
        plot(non_pml_x, non_pml_y, 'r--', 'LineWidth', 2);
        hold off;
    end
end

%normalize E_z, PML reflection, dx analyze, slit idx design

function [error_metrics] = calculate_error_metrics(E_z, E_z_analytical, PML)
    % Extract the non-PML region for fair comparison
    E_z_interior = E_z(:, PML+1:end-PML, PML+1:end-PML);
    E_z_analytical_interior = E_z_analytical(:, PML+1:end-PML, PML+1:end-PML);
    
    % Define a specific point to analyze (relative to interior grid)
    % Assuming we want point (150,150) in the global grid
    point_x = 160 - PML;
    point_y = 160 - PML;
    
    % Check if the point is within bounds of the interior grid
    interior_size = size(E_z_interior);
    if point_x > 0 && point_x <= interior_size(2) && point_y > 0 && point_y <= interior_size(3)
        valid_point = true;
    else
        valid_point = false;
        warning('Specified point (150,150) is outside the interior region. Using center point instead.');
        point_x = round(interior_size(2)/2);
        point_y = round(interior_size(3)/2);
    end
    
    % Calculate global error metrics at each time step
    for n = 1:size(E_z, 1)
        % L1 norm (mean absolute error)
        L1_error(n) = mean(abs(E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)), 'all');
        
        % L2 norm (root mean square error)
        L2_error(n) = sqrt(mean((E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)).^2, 'all'));
        
        % L∞ norm (maximum absolute error)
        Linf_error(n) = max(abs(E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)), [], 'all');
        
        % Calculate error at specific point (150,150)
        point_error(n) = abs(E_z_interior(n, point_x, point_y) - E_z_analytical_interior(n, point_x, point_y));
    end
    
    % Return all metrics
    error_metrics.L1 = L1_error;
    error_metrics.L2 = L2_error;
    error_metrics.Linf = Linf_error;
    error_metrics.point = point_error;
    
    % Calculate steady-state errors (use last ~80% of simulation)
    steady_idx = round(0.2*length(L1_error)):length(L1_error);
    error_metrics.steady_L1 = mean(L1_error(steady_idx));
    error_metrics.steady_L2 = mean(L2_error(steady_idx));
    error_metrics.steady_Linf = mean(Linf_error(steady_idx));
    error_metrics.steady_point = mean(point_error(steady_idx));
    
    % Create two subplots: global errors and point-specific error
    figure;
    
    % First subplot: Global error metrics
    subplot(2,1,1);
    plot(L1_error, 'b-', 'LineWidth', 2); hold on;
    plot(L2_error, 'r-', 'LineWidth', 2);
    plot(Linf_error, 'g-', 'LineWidth', 2);
    xline(steady_idx(1), 'k--', 'Steady State', 'LineWidth', 1.5);
    xlabel('Time Steps');
    ylabel('Error Magnitude');
    title('Global Error Metrics');
    legend('L1 (Mean Abs)', 'L2 (RMSE)', 'L∞ (Maximum)');
    grid on;
    
    % Second subplot: Point-specific error
    subplot(2,1,2);
    plot(point_error, 'c-', 'LineWidth', 2); hold on;
    xline(steady_idx(1), 'k--', 'Steady State', 'LineWidth', 1.5);
    xlabel('Time Steps');
    ylabel('Error Magnitude');
    if valid_point
        title(sprintf('Error at Point (150,150)'));
    else
        title(sprintf('Error at Center Point (grid limitations)'));
    end
    grid on;
    
    % Add text showing steady-state error values
    text(round(0.7*length(point_error)), max(point_error(steady_idx))*1.1, ...
        sprintf('Steady Point Error: %.2e', error_metrics.steady_point), 'Color', 'c');
    
    % Create separate figure for field value comparison at the point
    figure;
    E_z_point = squeeze(E_z_interior(:, point_x, point_y));
    E_z_analytical_point = squeeze(E_z_analytical_interior(:, point_x, point_y));
    
    plot(E_z_point, 'b-', 'LineWidth', 2); hold on;
    plot(E_z_analytical_point, 'r--', 'LineWidth', 2);
    xline(steady_idx(1), 'k--', 'Steady State', 'LineWidth', 1.5);
    xlabel('Time Steps');
    ylabel('E_z Amplitude');
    title('Field Value Comparison at Point (150,150)');
    legend('FDTD Solution', 'Analytical Solution');
    grid on;
    
    % Zoom in on steady-state region with another figure
    figure;
    plot(steady_idx, E_z_point(steady_idx), 'b-', 'LineWidth', 2); hold on;
    plot(steady_idx, E_z_analytical_point(steady_idx), 'r--', 'LineWidth', 2);
    xlabel('Time Steps');
    ylabel('E_z Amplitude');
    title('Field Values at Point (150,150) - Steady State Region');
    legend('FDTD Solution', 'Analytical Solution');
    grid on;
end

function analyze_pml_reflections(E_z, x, y, t, PML)
    % Choose time indices after waves have reached the PML
    start_idx = floor(length(t) * 0.3);
    end_idx = floor(length(t) * 0.8);
    
    % Define measurement lines near the PML
    measure_dist = 5; % cells away from PML
    
    % Initialize arrays to store maximum field values
    incident_max = zeros(1, length(t));
    reflected_max = zeros(1, length(t));
    
    % For each time step in our analysis window
    for n = start_idx:end_idx
        % Get field at current time
        field = squeeze(E_z(n, :, :));
        
        % Measure incident wave (near source, before PML reflection)
        center_x = round(length(x)/2);
        center_y = round(length(y)/2);
        incident_region = field(center_x-10:center_x+10, center_y-10:center_y+10);
        incident_max(n) = max(abs(incident_region(:)));
        
        % Measure potential reflections (just inside the PML boundary)
        reflection_region = field(PML+measure_dist:PML+2*measure_dist, PML+measure_dist:end-PML-measure_dist);
        reflected_max(n) = max(abs(reflection_region(:)));
    end
    
    % Calculate reflection coefficient (approximate)
    max_incident = max(incident_max);
    max_reflected = max(reflected_max(floor(length(t)*0.6):end)); % Look at later times for reflections
    reflection_coef = max_reflected / max_incident;
    
    % Plot results
    figure;
    semilogy(t, incident_max, 'b-', 'LineWidth', 2); hold on;
    semilogy(t, reflected_max, 'r-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Field Magnitude (log scale)');
    title(sprintf('PML Reflection Analysis (Coefficient ≈ %.2e)', reflection_coef));
    legend('Incident Wave', 'Reflected Wave');
    grid on;
    
    % Display assessment
    if reflection_coef < 1e-2
        disp('PML Performance: EXCELLENT (Reflection < 1%)');
    elseif reflection_coef < 5e-2
        disp('PML Performance: GOOD (Reflection < 5%)');
    elseif reflection_coef < 1e-1
        disp('PML Performance: ACCEPTABLE (Reflection < 10%)');
    else
        disp('PML Performance: POOR (Reflection > 10%)');
    end
end