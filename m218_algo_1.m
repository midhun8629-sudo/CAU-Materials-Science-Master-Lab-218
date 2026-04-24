clear all; 
clc;

% 1. Define Physical and Numerical Parameters
l = 1.0;          
D = 1.0;          
N = 5000;           
dx = l / N;       

% Choose dt to satisfy the explicit stability criterion: dt <= dx^2 / (2*D)
dt = (dx^2) / (2.1 * D); 
M = 1000;         
T = M * dt;       

% 2. Initialize Arrays (Task 4.1 I a-c)
x = linspace(0, l, N+1); 
t = linspace(0, T, M+1); 
c = zeros(N+1, M+1);     

% 3. Apply Initial Condition at t=0, first column of c. select all rows
% within the first column
c(:, 1) =  exp(-50 * (x/l - 0.5).^2); 

% 4. Apply Boundary Conditions at t=0
% c(1, 1) = c(2,1);    % Derichlet prescribe conc      
%c(N+1, 1) = c(N, 1); % Neumann zero heat flux  

% 5. Explicit Time Integration Loop 
for i = 1:M
    % Apply Dirichlet boundary condition at the new time step. first
    % first row of c 
    %c(1, i+1) = c(2,i+1);
    % Compute interior spatial nodes
    for n = 2:N
        c(n, i+1) = c(n, i) + (D * dt / dx^2) * (c(n+1, i) - 2*c(n, i) + c(n-1, i));
    end
    
    % Apply Neumann boundary condition at the new time step
    c(N+1, i+1) = c(N, i+1); 
    c(1, i+1) = c(2,i+1);

end

% 6. Visualization
figure;
plotFreq = 50; 
pauseTime = 0.05;

for i = 1:plotFreq:M+1
    plot(x, c(:, i), 'r', 'LineWidth', 2);
    axis([0 1 0 1]);
    xlabel('Location x');
    ylabel('Concentration c');
    title(sprintf('Time step %d', i-1));
    grid on;
    pause(pauseTime);
end
% 7. Surface Plot of Concentration over Space and Time (Task 4.1 II b)
figure; 
[T_grid, X_grid] = meshgrid(t, x); 
surf(T_grid, X_grid, c, 'EdgeColor', 'none'); 

colormap jet; 
colorbar; 

xlabel('Time t');
ylabel('Location x');
zlabel('Concentration c');
title('Spatiotemporal Evolution of Concentration');
view(-37.5, 30); 
grid on;