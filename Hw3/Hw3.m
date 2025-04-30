% Constants
epsilon0 = 8.854187817e-12; % Vacuum permittivity (F/m)
L = 0.01; % Plate size (m)

% List of N values to test
N_list = [10, 20, 30, 50, 100];

for idx = 1:length(N_list)
    N = N_list(idx);
    fprintf('\n--- N = %d x %d ---\n', N, N);

    a = L / N; % Cell size
    A = a^2;   % Cell area

    % Generate center points
    centers = zeros(N*N, 2);
    count = 1;
    for i = 0:N-1
        for j = 0:N-1
            centers(count,1) = (i+0.5)*a - L/2;
            centers(count,2) = (j+0.5)*a - L/2;
            count = count + 1;
        end
    end

    % Initialize Z matrix
    Z = zeros(N*N);

    % Fill Z matrix
    % Initialize Z matrix for saving results
    Z_all = cell(length(N_list), 1);

    for i = 1:N*N
        xi = centers(i,1);
        yi = centers(i,2);
        for j = 1:N*N
            xj = centers(j,1);
            yj = centers(j,2);
            if i == j
                % Self-term integration
                integrand = @(x,y) 1 ./ sqrt(x.^2 + y.^2);
                result = integral2(integrand, -a/2, a/2, @(x) -a/2 * ones(size(x)), @(x) a/2 * ones(size(x)));
                Z(i,j) = result / (4*pi*epsilon0);
            else
                Rij = sqrt((xi - xj)^2 + (yi - yj)^2);
                Z(i,j) = A / (4*pi*epsilon0*Rij);
            end
        end
    end



    % Solve for surface charge density
    V = ones(N*N, 1); % Voltage = 1 V everywhere
    sigma = Z\V;

    % Average self-term
    avg_Zii = mean(diag(Z));
    fprintf('Average Z_ii: %.3e Ohm\n', avg_Zii);

    % Total charge and capacitance
    Q = sum(sigma)*A;
    C = Q / 1; % Since V = 1V
    fprintf('Capacitance: %.3e F\n', C);

    % Plot charge distribution
    sigma_2D = reshape(sigma, [N, N]);
    figure;
    imagesc(linspace(-L/2,L/2,N), linspace(-L/2,L/2,N), sigma_2D);
    colorbar;
    c = colorbar;
    c.Label.String = 'Charge Density (C/m^2)'; % Label for colorbar
    title(sprintf('Charge Distribution (N=%d)', N));
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;
end
