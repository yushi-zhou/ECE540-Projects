% A Matlab program using FDTD to calculate the radiation of a infinite current line, with Yee's grid.
% Author: Yushi Zhou  
function infinite_current()
    %define physical constants
    c = 3e8; %speed of light
    mu0 = 4*pi*1e-7; %permeability of free space
    epsilon0 = 8.85e-12; %permittivity of free space
    freq = 1e9;
    lambda = c/freq;
    x_source = 1;
    y_source = 1;

    % Define the grid
    dx = 0.005;
    dy = dx;
    dt = dx/5e8; %  dt<=\frac{dx}{c\sqrt{2}}

    %define range
    x = 0:dx:2;
    y = 0:dy:2;
    %define time range
    t = 0:dt:15/freq;

    %define permitivity and conductivity
    epsilon = ones(length(x),length(y))*epsilon0;
    sigma = zeros(length(x),length(y));

    %define a sin current source, period = 100*dt, (1.5,1)
    omega = 2*pi*freq;
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
    
    % plot
    visualize_field(E_z, x, y, t, dt);
 
    save('infinite_current.mat','E_z','H_x','H_y','J_z','x','y','t');

end

function E_z_analytical = analytical_solution(r, t, omega, mu0, epsilon0)
    %Helmholtz solution
    k = omega * sqrt(mu0 * epsilon0);
    E_z_analytical = (1i/4)*besselh(0, 2, k*r - omega*t); % Hankel function of second kind
    E_z_analytical = real(E_z_analytical); % Take real part
end
    


function visualize_field(E_z, x, y, t, dt)
    figure;
    hFig = gcf;
    hAx = axes('Parent', hFig);
    hImg = imagesc(x, y, squeeze(E_z(1, :, :))', 'Parent', hAx);
    colorbar;
    colormap(winter);
    %clim([-0.8, 0.8]); % Set color scale limits for better contrast
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