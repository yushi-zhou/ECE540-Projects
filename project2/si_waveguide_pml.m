clear all; close all; clc;

a = 0.033;    
b = 0.014;   
L = a/3;      % air+pml thickness

lambda0 = 1.55e-2;  % [m]
k0 = 2*pi/lambda0;  
alpha = k0;         % PML absorption parameter (tunable)

R1 = [3; 4; 0;    a;  a;  0;  0;  0;  b;  b];
R2 = [3; 4; -L; a+L; a+L; -L; -L; -L; b+L; b+L];

gd = [R1, R2];
ns = char('R1','R2')';
sf = 'R2';
dl = decsg(gd, sf, ns);

model = createpde();
geometryFromEdges(model, dl);

hmax = a/100;
generateMesh(model, 'Hmax', hmax, 'GeometricOrder', 'linear');

p = model.Mesh.Nodes;
t = model.Mesh.Elements;
numNodes = size(p, 2);    
numElements = size(t, 2);

figure;
pdemesh(model);
axis equal;
title('Mesh for Si Waveguide with Air & PML Regions');
xlabel('x (m)'); ylabel('y (m)');

hold on;
waveguide_x = [0, a, a, 0, 0];
waveguide_y = [0, 0, b, b, 0];
plot(waveguide_x, waveguide_y, 'r-', 'LineWidth', 2);

pml_thickness = L*0.3;
air_boundary_x = [-L+pml_thickness, a+L-pml_thickness, a+L-pml_thickness, -L+pml_thickness, -L+pml_thickness];
air_boundary_y = [-L+pml_thickness, -L+pml_thickness, b+L-pml_thickness, b+L-pml_thickness, -L+pml_thickness];
plot(air_boundary_x, air_boundary_y, 'b--', 'LineWidth', 2);

pml_outer_x = [-L, a+L, a+L, -L, -L];
pml_outer_y = [-L, -L, b+L, b+L, -L];
plot(pml_outer_x, pml_outer_y, 'k-', 'LineWidth', 2);

legend('Mesh', 'Waveguide', 'Air-PML Boundary', 'Simulation Boundary', 'Location', 'northeastoutside');
hold off;

%Si and air permitivity
epsilon_Si = 12.04;
epsilon_air = 1;

elemEpsilon = zeros(1, numElements);
for e = 1:numElements
    nodes = t(:, e);
    x_coords = p(1, nodes);
    y_coords = p(2, nodes);
    xc = mean(x_coords);
    yc = mean(y_coords);
    if (xc >= 0 && xc <= a && yc >= 0 && yc <= b)
        elemEpsilon(e) = epsilon_Si;
    else
        elemEpsilon(e) = epsilon_air;
    end
end

%Matrix Assembly
K = sparse(numNodes, numNodes);
M = sparse(numNodes, numNodes);

for e = 1:numElements
    nodes = t(:, e);
    p_e = p(:, nodes);
    x = p_e(1, :);
    y = p_e(2, :);
    
    A_e = abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))/2;
    
    b_i = [y(2)-y(3); x(3)-x(2)];
    b_j = [y(3)-y(1); x(1)-x(3)];
    b_k = [y(1)-y(2); x(2)-x(1)];
    
    K_e = (1/(4*A_e)) * [ ...
        dot(b_i, b_i), dot(b_i, b_j), dot(b_i, b_k); ...
        dot(b_j, b_i), dot(b_j, b_j), dot(b_j, b_k); ...
        dot(b_k, b_i), dot(b_k, b_j), dot(b_k, b_k) ];
    
    epsilon_local = elemEpsilon(e);
    M_e = epsilon_local*(A_e/12)*[2,1,1; 1,2,1; 1,1,2];
    
    for i = 1:3
        for j = 1:3
            K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + K_e(i,j);
            M(nodes(i), nodes(j)) = M(nodes(i), nodes(j)) + M_e(i,j);
        end
    end
end

%% PML
pml_thickness = L *0.3;  % PML thickness
sigma_max = 10;          % Maximum PML strength
pml_order = 3;           % Polynomial order for PML profile

sigma_x = zeros(1, numElements);
sigma_y = zeros(1, numElements);

for e = 1:numElements
    nodes = t(:, e);
    x_coords = p(1, nodes);
    y_coords = p(2, nodes);
    xc = mean(x_coords);
    yc = mean(y_coords);
    
    % PML in x-direction (left and right boundaries)
    if xc < -L + pml_thickness  % Left PML
        d = -L + pml_thickness - xc;  % Distance into PML
        sigma_x(e) = sigma_max * (d/pml_thickness)^pml_order;
    elseif xc > a + L - pml_thickness  % Right PML
        d = xc - (a + L - pml_thickness);  % Distance into PML
        sigma_x(e) = sigma_max * (d/pml_thickness)^pml_order;
    end
    
    % PML in y-direction (bottom and top boundaries)
    if yc < -L + pml_thickness  % Bottom PML
        d = -L + pml_thickness - yc;  % Distance into PML
        sigma_y(e) = sigma_max * (d/pml_thickness)^pml_order;
    elseif yc > b + L - pml_thickness  % Top PML
        d = yc - (b + L - pml_thickness);  % Distance into PML
        sigma_y(e) = sigma_max * (d/pml_thickness)^pml_order;
    end
end

% Rebuild matrices with PML
K = sparse(numNodes, numNodes);
M = sparse(numNodes, numNodes);

for e = 1:numElements
    nodes = t(:, e);
    p_e = p(:, nodes);
    x = p_e(1, :);
    y = p_e(2, :);
    
    A_e = abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))/2;
    
    b_i = [y(2)-y(3); x(3)-x(2)];
    b_j = [y(3)-y(1); x(1)-x(3)];
    b_k = [y(1)-y(2); x(2)-x(1)];
    
    % PML stretching factors for this element
    sx = 1 + 1i*sigma_x(e)/k0;  % Complex coordinate stretch in x
    sy = 1 + 1i*sigma_y(e)/k0;  % Complex coordinate stretch in y
    
    %local K
    K_e = (1/(4*A_e)) * [ ...
        (b_i(1)^2/sx + b_i(2)^2/sy), (b_i(1)*b_j(1)/sx + b_i(2)*b_j(2)/sy), (b_i(1)*b_k(1)/sx + b_i(2)*b_k(2)/sy); ...
        (b_j(1)*b_i(1)/sx + b_j(2)*b_i(2)/sy), (b_j(1)^2/sx + b_j(2)^2/sy), (b_j(1)*b_k(1)/sx + b_j(2)*b_k(2)/sy); ...
        (b_k(1)*b_i(1)/sx + b_k(2)*b_i(2)/sy), (b_k(1)*b_j(1)/sx + b_k(2)*b_j(2)/sy), (b_k(1)^2/sx + b_k(2)^2/sy)];
    
    % Mass matrix with material weight and PML stretch 
    epsilon_local = elemEpsilon(e);
    det_s = sx * sy;  % mapping coefficiency
    M_e = epsilon_local * (A_e/12) * det_s * [2,1,1; 1,2,1; 1,1,2];
    
    for i = 1:3
        for j = 1:3
            K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + K_e(i,j);
            M(nodes(i), nodes(j)) = M(nodes(i), nodes(j)) + M_e(i,j);
        end
    end
end

A_mat = K;

%solve
numModes = 10;
[V, D] = eigs(A_mat, M, numModes, 'smallestabs');
beta_squared = diag(D);

[beta_squared, idx] = sort(beta_squared);
V = V(:, idx);

beta_squared = beta_squared(2:10);
V = V(:, 2:10);
numModes = 9;

beta = sqrt(beta_squared);
neff = beta ./ k0;

fprintf('Si Waveguide Simulation Results:\n');
fprintf('-------------------------------------------\n');
fprintf('Mode  |    β (rad/m)    |   Effective Index\n');
fprintf('-------------------------------------------\n');
for i = 1:numModes
    fprintf(' %d    |   %10.2f   |    %8.4f\n', i, beta(i), neff(i));
end

%plot
figure;
nRows = ceil(sqrt(numModes));
nCols = ceil(numModes/nRows);

air_boundary_x = [-L+pml_thickness, a+L-pml_thickness, a+L-pml_thickness, -L+pml_thickness, -L+pml_thickness];
air_boundary_y = [-L+pml_thickness, -L+pml_thickness, b+L-pml_thickness, b+L-pml_thickness, -L+pml_thickness];

waveguide_x = [0, a, a, 0, 0];
waveguide_y = [0, 0, b, b, 0];

pml_outer_x = [-L, a+L, a+L, -L, -L];
pml_outer_y = [-L, -L, b+L, b+L, -L];

for i = 1:numModes
    subplot(nRows, nCols, i);
    field = real(V(:, i));
    trisurf(t(1:3,:)', p(1,:), p(2,:), field, 'EdgeColor','none');
    view(2);
    
    hold on;
    plot3(air_boundary_x, air_boundary_y, max(field)*ones(size(air_boundary_x)), 'r-', 'LineWidth', 1.5);
    plot3(pml_outer_x, pml_outer_y, max(field)*ones(size(pml_outer_x)), 'k--', 'LineWidth', 1);
    plot3(waveguide_x, waveguide_y, max(field)*ones(size(waveguide_x)), 'w-', 'LineWidth', 1.5);
    hold off;
    
    axis equal tight;
    title(sprintf('Mode %d: n_eff=%.4f', i, neff(i)));
    xlabel('x (m)'); ylabel('y (m)');
    colorbar;
    colormap(jet);
    
    if i == 1
        legend('', 'Air-PML boundary', 'Simulation boundary', 'Waveguide', 'Location', 'northeastoutside');
    end
end

sgtitle('Si Waveguide Modes - FEM Simulation with PML');
