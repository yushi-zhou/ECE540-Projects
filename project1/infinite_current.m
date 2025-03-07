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
    sigma_max = -log(0.001)*5/(2*377*PML*dx);
    for i = 1:PML
        for j = 1:length(y)
            sigma(i, j) = sigma_max * ((PML - i) / PML)^4;
            sigma(length(x) - i + 1, j) = sigma_max * ((PML - i) / PML)^4;
        end
    end
    for i = 1:length(x)
        for j = 1:PML
            sigma(i, j) = sigma_max * ((PML - j) / PML)^3;
            sigma(i, length(y) - j + 1) = sigma_max * ((PML - j) / PML)^4;
        end
    end

    % Define PML at the corners
    for i = 1:PML
        for j = 1:PML
            sigma_x = sigma_max * ((PML - i + 0.5) / PML)^4;
            sigma_y = sigma_max * ((PML - j + 0.5) / PML)^4;
            
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
    %error_metrics = calculate_error_metrics(E_z, E_z_analytical, PML, dx, dy);
    
    
    % plot
    visualize_field(E_z, x, y, t, PML);
 
    save('infinite_current.mat','E_z','H_x','H_y','J_z','x','y','t');

    % Calculate PML reflection
    calculate_pml_reflection(E_z, x, y, t, PML);
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
    clim([-1000, 1000]); % Set color scale limits for better contrast
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
        title(hAx, sprintf('E_z at t = %.2e s', t(n)));
        
        % Ensure PML boundary remains visible
        hold on;
        plot(non_pml_x, non_pml_y, 'r--', 'LineWidth', 2);
        hold off;
    end
end

%normalize E_z, PML reflection, dx analyze, slit idx design

function [error_metrics] = calculate_error_metrics(E_z, E_z_analytical, PML, dx, dy)
    % Extract the non-PML region for fair comparison
    E_z_interior = E_z(:, PML+1:end-PML, PML+1:end-PML);
    E_z_analytical_interior = E_z_analytical(:, PML+1:end-PML, PML+1:end-PML);
    
    % Define a point to analyze
    x_observe = 0.5;
    y_observe = 0.5;
    point_x = round(x_observe / (dx));
    point_y = round(y_observe / (dy));
    
    % Check if the point is within bounds of the interior grid
    interior_size = size(E_z_interior);
    
    % Calculate global error metrics at each time step
    for n = 1:size(E_z, 1)
        % Calculate error at specific point 
        point_error(n) = abs((E_z_interior(n, point_x, point_y) - E_z_analytical_interior(n, point_x, point_y)));
    end
    
    % Return metrics
    error_metrics.point = point_error;
    
    % Calculate steady-state errors 
    steady_idx = round(0.8*length(point_error)):length(point_error);
    error_metrics.steady_point = mean(point_error(steady_idx));
    
    % Create figure: point-specific error
    figure;
    plot(point_error, 'r-', 'LineWidth', 2); hold on;
    %xline(steady_idx(1), 'k--', 'Steady State', 'LineWidth', 1.5);
    xlabel('Time Steps');
    ylabel('Error');
    title(sprintf('Error at Point (0.5,0.5)'));
    
    grid on;
    
    % Add text showing steady-state error values
    text(round(0.7*length(point_error)), max(point_error(steady_idx))*1.1, ...
        sprintf('Steady Point Error: %.2e', error_metrics.steady_point), 'Color', 'r');
    
    % Create separate figure for field value comparison at the point
    figure;
    E_z_point = squeeze(E_z_interior(:, point_x, point_y));
    E_z_analytical_point = squeeze(E_z_analytical_interior(:, point_x, point_y));
    
    plot(E_z_point, 'b-', 'LineWidth', 2); hold on;
    plot(E_z_analytical_point, 'r--', 'LineWidth', 2);
    %xline(steady_idx(1), 'k--', 'Steady State', 'LineWidth', 1.5);
    xlabel('Time Steps');
    ylabel('E_z Amplitude');
    title('Field Value Comparison at Point (0.5,0.5)');
    legend('FDTD Solution', 'Analytical Solution');
    grid on;
end

function calculate_pml_reflection(E_z, x, y, t, PML)
    % Calculate reflection from PML by comparing field values near boundary
    % Place observation points: one just before PML, one inside computational domain
    obs_inside = PML + 10; % 10 cells inside computational domain
    obs_near_pml = PML + 3;  % 3 cells from PML boundary
    mid_y = round(length(y)/2);
    
    % Extract time-domain signals at these points
    signal_inside = squeeze(E_z(:, obs_inside, mid_y));
    signal_near_pml = squeeze(E_z(:, obs_near_pml, mid_y));
    
    % Wait for steady state (skip initial transients)
    steady_start = round(0.7*length(t));
    
    % Calculate reflection coefficient from ratio of amplitudes
    % For perfect PML, waves at the boundary should pass through without reflection
    % So the amplitude ratio should approach 1 in steady state
    max_inside = max(abs(signal_inside(steady_start:end)));
    max_near_pml = max(abs(signal_near_pml(steady_start:end)));
    
    % Simple reflection estimate
    reflection_coeff = abs(max_near_pml/max_inside - 1);
    
    % Plot the time-domain signals
    figure;
    plot(t, signal_inside, 'b-', 'LineWidth', 2); hold on;
    plot(t, signal_near_pml, 'r--', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('E_z Amplitude');
    title('PML Reflection Analysis');
    legend('10 cells from PML', '3 cells from PML');
    grid on;
    
    % Add text showing reflection coefficient
    text(t(round(length(t)*0.7)), max(abs(signal_inside))*0.8, ...
        sprintf('Est. Reflection: %.2e', reflection_coeff), 'FontSize', 12);
    
    fprintf('Estimated PML reflection coefficient: %.2e\n', reflection_coeff);
    
end

