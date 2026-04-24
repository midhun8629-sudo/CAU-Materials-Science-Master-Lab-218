clear all; clc; close all;

%% ========================================================================
%  SCENARIO SWITCHBOARD
%  Uncomment ONLY ONE line below to run the specific experiment.
% ========================================================================

% --- Explicit Experiments (Section 4.1) ---
% scenario = 1; % Explicit: Standard Gaussian (Stable Animation) 
% scenario = 2; % Explicit: Instability (Violate dt limit)
% scenario = 5; % Explicit: Block IC (Static Subplots: t=0, 300, 1000)
% scenario = 6; % Explicit: Double Dirichlet Boundary (Animation)

% --- Implicit Experiments (Section 4.2) ---
% scenario = 3;  % Implicit: Standard Gaussian (Stable Animation) 
% scenario = 4;  % Implicit: Unconditional Stability (Massive dt)
scenario = 9;  % Implicit: Block IC (Static Subplots: t=0, 300, 1000)
% scenario = 10; % Implicit: Double Dirichlet (Static Subplots: t=0, 15, 30)

% --- Lab Report Exact Figure Generation ---
% scenario = 7; % Explicit: Fig. 4 Static Subplots (t=0, 300, 1000)
% scenario = 8; % Implicit: Fig. 6 Static Subplots (t=0, 15, 30)

%% 1. Define Base Physical Parameters
l = 1.0;          
D = 1.0;          
N = 50;           
dx = l / N;       

%% 2. Scenario Configuration Logic
if ismember(scenario, [1, 5, 6, 7])
    solver = 'explicit';
    dt = (dx^2) / (2.1 * D); 
    M = 1000;
elseif scenario == 2
    solver = 'explicit';
    dt = (dx^2) / (1.9 * D); % UNSTABLE
    M = 200; 
elseif ismember(scenario, [3, 9, 10])
    solver = 'implicit';
    dt = (dx^2) / (2.1 * D); % Use standard dt for fair comparison
    M = 1000;
elseif scenario == 4
    solver = 'implicit';
    dt = 0.05; % MASSIVE dt 
    M = 50;                  
elseif scenario == 8
    solver = 'implicit';
    dt = 0.01; 
    M = 30;                  
end

T = M * dt;       

%% 3. Initialize Arrays
x = linspace(0, l, N+1)'; 
t = linspace(0, T, M+1); 
c = zeros(N+1, M+1);     

%% 4. Apply Initial Condition at t=0
if ismember(scenario, [5, 9])
    % Modified IC: Block / Top-Hat Distribution
    c(:, 1) = 0; 
    c(round(N/3):round(2*N/3), 1) = 1.0; 
else
    % Standard IC: Gaussian Curve
    c(:, 1) = exp(-50 * (x/l - 0.5).^2); 
end

%% 5. Time Integration
if strcmp(solver, 'explicit')
    % --- EXPLICIT SOLVER ---
    c(1, 1) = 0; 
    if ismember(scenario, [6, 10]) % Applying to 6 just in case
        c(N+1, 1) = 0.8; 
    else
        c(N+1, 1) = c(N, 1); 
    end
    
    for i = 1:M
        c(1, i+1) = 0; 
        for n = 2:N
            c(n, i+1) = c(n, i) + (D * dt / dx^2) * (c(n+1, i) - 2*c(n, i) + c(n-1, i));
        end
        if ismember(scenario, [6, 10])
            c(N+1, i+1) = 0.8; 
        else
            c(N+1, i+1) = c(N, i+1); 
        end
    end
    
elseif strcmp(solver, 'implicit')
    % --- IMPLICIT SOLVER ---
    c(1, 1) = 0; 
    if scenario == 10
        c(N+1, 1) = 0.8; % Initial state for Right Dirichlet
    else
        c(N+1, 1) = c(N, 1);
    end
    
    a_val = -D * dt / (dx^2);
    b_val = 1 - 2 * a_val;   
    
    % Adjust diagonals based on boundary condition
    if scenario == 10
        A = [a_val * ones(1, N-1), 0];   % 0 at the end prevents neighbor coupling 
    else
        A = [a_val * ones(1, N-1), -1];  % -1 at the end enforces Neumann derivative
    end
    
    B = [1, b_val * ones(1, N-1), 1];    
    C = [0, a_val * ones(1, N-1)];       
    
    for i = 1:M
        R = c(:, i); 
        R(1) = 0;       
        if scenario == 10
            R(end) = 0.8; % Force Right Dirichlet value in RHS vector
        else
            R(end) = 0;   % Force Right Neumann zero-flux in RHS vector
        end
        
        [P, G] = triDiagTrafo(R, A, B, C);
        c(:, i+1) = triDiagSolve(c(:, i), P, G);
    end
end

%% 6. Visualization: 2D Plots
% Scenarios that generate animated exploratory plots:
if ismember(scenario, [1, 2, 3, 4, 6])
    figure(1);
    plotFreq = max(1, round(M/20)); 
    pauseTime = 0.05;

    for i = 1:plotFreq:M+1
        plot(x, c(:, i), 'r-', 'LineWidth', 2);
        axis([0 l -0.1 1.2]); 
        xlabel('Location x');
        ylabel('Concentration c');
        title(sprintf('Scenario %d | Time step %d | Solver: %s', scenario, i-1, solver));
        grid on;
        pause(pauseTime);
    end

% Scenarios that generate static 1x3 subplots for the lab report:
elseif ismember(scenario, [5, 7, 8, 9, 10])
    figure('Name', sprintf('Static Subplots - Scenario %d', scenario), 'Position', [100, 100, 1200, 400]);
    
    if ismember(scenario, [5, 7]) % Explicit requires higher time steps
        t_steps = [1, 301, 1001]; 
        titles = {'time step 0', 'time step 300', 'time step 1000'};
    elseif ismember(scenario, [8, 9, 10]) % Implicit comparisons
        t_steps = [1, 301, 1001];    
        titles = {'time step 0', 'time step 300', 'time step 1000'};
    end
    
    for k = 1:3
        subplot(1, 3, k);
        plot(x, c(:, t_steps(k)), 'r-', 'LineWidth', 2);
        axis([0 l -0.1 1.2]); % 1.2 height accommodates the Block IC perfectly
        xlabel('location x');
        ylabel('concentration c');
        title(titles{k});
        grid on;
    end
end

%% 7. Visualization: 3D Surface Plot
if ismember(scenario, [1, 2, 3, 4, 5, 6, 9, 10])
    figure(2); 
    [T_grid, X_grid] = meshgrid(t, x); 
    surf(T_grid, X_grid, c, 'EdgeColor', 'none'); 
    colormap jet; 
    colorbar; 
    xlabel('Time t');
    ylabel('Location x');
    zlabel('Concentration c');
    title(sprintf('Spatiotemporal Evolution | Scenario %d', scenario));
    view(-37.5, 30); 
    grid on;
end

%% ==========================================
%  Local Functions (Thomas Algorithm)
% ==========================================
function [P, G] = triDiagTrafo(R, A, B, C)
    N_plus_1 = length(B);
    P = zeros(N_plus_1, 1); 
    G = zeros(N_plus_1 - 1, 1);
    P(1) = R(1) / B(1);
    G(1) = C(1) / B(1);
    for i = 2:N_plus_1
        denom = B(i) - A(i-1) * G(i-1); 
        P(i) = (R(i) - A(i-1) * P(i-1)) / denom;
        if i < N_plus_1
            G(i) = C(i) / denom;
        end
    end
end

function c_new = triDiagSolve(~, P, G)
    N_plus_1 = length(P);
    c_new = zeros(N_plus_1, 1);
    c_new(N_plus_1) = P(N_plus_1);
    for i = (N_plus_1 - 1):-1:1
        c_new(i) = P(i) - G(i) * c_new(i+1);
    end
end