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
        J_z(n,x_source/dx,y_source/dy) = sin(omega*n*dt);
    end

    %define perfect conductor sheet at x=x_conductor
    idx_conductor_x = 50;
    idx_slit_y = 0.3*length(y);

    %define PML 
    PML = (lambda)/dx;
    sigma_max = 0.1;
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
            sigma(i, j) = sigma_max * (((PML - i) / PML)^3 + ((PML - j) / PML)^3);
            sigma(i, length(y) - j + 1) = sigma_max * (((PML - i) / PML)^3 + ((PML - j) / PML)^3);
            sigma(length(x) - i + 1, j) = sigma_max * (((PML - i) / PML)^3 + ((PML - j) / PML)^3);
            sigma(length(x) - i + 1, length(y) - j + 1) = sigma_max * (((PML - i) / PML)^3 + ((PML - j) / PML)^3);
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
                E_z_analytical(n, i, j) = analytical_solution(r, n*dt, omega, mu0, epsilon0);
            end
        end
    end
    error = abs(E_z - E_z_analytical); %error at the beginning is large because system is not in steady yet
    
    % Add this after calculating 'error'
    error_metrics = calculate_error_metrics(E_z, E_z_analytical, PML);
    
    % plot
    visualize_field(E_z, x, y, t, dt);
    visualize_field(E_z_analytical, x, y, t, dt);
    visualize_field(error, x, y, t, dt);
 
    save('infinite_current.mat','E_z','H_x','H_y','J_z','x','y','t');

end

function E_z_analytical = analytical_solution(r, t, omega, mu0, epsilon0)
    %Helmholtz solution
    k = omega * sqrt(mu0 * epsilon0);
    E_z_analytical = (1/4)*besselh(0, 1, -k*r )*exp(1i*omega*t); % Hankel function of second kind
    E_z_analytical = imag(E_z_analytical); 
end
    


function visualize_field(E_z, x, y, t, dt)
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

    % Create a slid bar
    hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(t), 'Value', 1, ...
                        'Units', 'normalized', 'Position', [0.2 0.01 0.6 0.05], ...
                        'Callback', @(src, event) updateImage());

    % update the image with the bar
    addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateImage());
    function updateImage()
        n = round(get(hSlider, 'Value'));
        set(hImg, 'CData', squeeze(E_z(n, :, :))');
        title(hAx, sprintf('E_z at t = %.2f', t(n)/dt));
    end
end

%normalize E_z, PML reflection, dx analyze, slit idx design

function [error_metrics] = calculate_error_metrics(E_z, E_z_analytical, PML)
    % Extract the non-PML region for fair comparison
    E_z_interior = E_z(:, PML+1:end-PML, PML+1:end-PML);
    E_z_analytical_interior = E_z_analytical(:, PML+1:end-PML, PML+1:end-PML);
    
    % Calculate error at each time step
    for n = 1:size(E_z, 1)
        % L1 norm (mean absolute error)
        L1_error(n) = mean(abs(E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)), 'all');
        
        % L2 norm (root mean square error)
        L2_error(n) = sqrt(mean((E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)).^2, 'all'));
        
        % L∞ norm (maximum absolute error)
        Linf_error(n) = max(abs(E_z_interior(n,:,:) - E_z_analytical_interior(n,:,:)), [], 'all');
    end
    
    % Return all metrics
    error_metrics.L1 = L1_error;
    error_metrics.L2 = L2_error;
    error_metrics.Linf = Linf_error;
    
    % Calculate steady-state errors (use last ~80% of simulation)
    steady_idx = round(0.2*length(L1_error)):length(L1_error);
    error_metrics.steady_L1 = mean(L1_error(steady_idx));
    error_metrics.steady_L2 = mean(L2_error(steady_idx));
    error_metrics.steady_Linf = mean(Linf_error(steady_idx));
    
    % Plot error metrics
    figure;
    time = 1:length(L1_error);
    
    % Create plot with all three error metrics
    plot(time, L1_error, 'b-', 'LineWidth', 2); hold on;
    plot(time, L2_error, 'r-', 'LineWidth', 2);
    plot(time, Linf_error, 'g-', 'LineWidth', 2);
    
    % Highlight steady-state region
    xline(steady_idx(1), 'k--', 'Steady State Region', 'LineWidth', 1.5);
    
    % Add labels and legend
    xlabel('Time Steps');
    ylabel('Error Magnitude');
    title('Error Evolution Over Simulation Time');
    legend('L1 (Mean Absolute Error)', 'L2 (RMSE)', 'L∞ (Maximum Error)');
    grid on;
    
    % Add text showing steady-state error values
    text_x = round(0.7*length(L1_error));
    text_y_offset = 0.05*(max(Linf_error) - min(L1_error));
    text(text_x, max(L1_error(steady_idx))+text_y_offset, ...
        sprintf('Steady L1: %.2e', error_metrics.steady_L1), 'Color', 'blue');
    text(text_x, max(L2_error(steady_idx))+2*text_y_offset, ...
        sprintf('Steady L2: %.2e', error_metrics.steady_L2), 'Color', 'red');
    text(text_x, max(Linf_error(steady_idx))+3*text_y_offset, ...
        sprintf('Steady L∞: %.2e', error_metrics.steady_Linf), 'Color', 'green');
end