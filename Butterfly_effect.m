% Lorenz Attractor Simulation
% Parameters for the Lorenz attractor
sigma = 10;
rho = 28;
beta = 8/3;

% Time settings
dt = 0.01;            % Time step
t_end = 20;          % End time
t = 0:dt:t_end;      % Time vector

% Initial conditions
initial_conditions = [1, 1, 1;   % Initial state for the first trajectory
                      1.01, 1, 1; % Slightly perturbed initial state
                      1, 1.01, 1; % Another perturbation
                      1, 1, 1.01; % Another perturbation
                      1.01, 1.01, 1; % Combined perturbation
                      1.01, 1, 1.01; % Combined perturbation
                      1, 1.01, 1.01; % Combined perturbation
                      1.01, 1.01, 1.01]; % Combined perturbation

% Prepare figure
figure;
hold on;

% Iterate over initial conditions
for i = 1:size(initial_conditions, 1)
    % Initialize variables for the trajectory
    xyz = zeros(length(t), 3);
    xyz(1, :) = initial_conditions(i, :);
    
    % Calculate the trajectory using the Lorenz equations
    for j = 1:length(t)-1
        x = xyz(j, 1);
        y = xyz(j, 2);
        z = xyz(j, 3);
        
        % Lorenz equations
        dx = sigma * (y - x);
        dy = x * (rho - z) - y;
        dz = x * y - beta * z;
        
        % Update the next state
        xyz(j + 1, :) = xyz(j, :) + dt * [dx, dy, dz];
    end
    
    % Plot the trajectory
    plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'DisplayName', ['Initial Conditions ' num2str(i)]);
end

% Final touches on the plot
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Lorenz Attractor - Butterfly Effect');
grid on;
view(30, 10);  % Set the viewing angle

% Adjust legend position
lgd = legend('show');
lgd.Location = 'northeastoutside';  % Move legend outside the plot
lgd.Position(1) = lgd.Position(1) + 0.1; % Shift the legend box to the right
lgd.Position(2) = lgd.Position(2) + 0.1; % Shift the legend box up (if needed)

hold off;
