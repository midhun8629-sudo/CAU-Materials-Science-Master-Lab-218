clear all; 
clc;
% Implicit Scheme
% 1. Physical and Numerical Parameters
l = 1.0;          
D = 1.0;          
N = 50;           
dx = l / N;       

% Implicit solvers are unconditionally stable, allowing massive time steps
dt = 0.001; 
M = 100;         
T = M * dt;      

% 2. Initialization & Boundary Conditions
x = linspace(0, l, N+1)'; % Column vector
c = zeros(N+1, 1);        % Overwriting array to save memory (COLUMN VECTOR)

% Initial condition at t=0
c = exp(-50 * (x/l - 0.5).^2);

% Boundary conditions for the initial state
c(1) = 0; % derichlet         
c(N+1) = c(N);    %neumann

% 3. Construct Diagonals A, B, and C
a_val = -D * dt / (dx^2);
b_val = 1 - 2 * a_val;   

A = [a_val * ones(1, N-1), -1];      % Subdiagonal
B = [1, b_val * ones(1, N-1), 1];    % Main diagonal
C = [0, a_val * ones(1, N-1)];       % Superdiagonal

% 4. Implicit Time Integration Loop
% figure; hold on;
for i = 2:M+1
    % Set right-hand side to current concentration
    R = c; 
    % Enforce the mathematical boundary conditions on the RHS vector
    R(1) = 0;       % Dirichlet
    R(end) = 0;     % Neumann
    
    % Forward Sweep (Gaussian elimination)
    [P, G] = triDiagTrafo(R, A, B, C);
    
    % Backward Substitution to find new concentration
    c = triDiagSolve(c, P, G);
    
    % Visualization
    plot(x, c, 'b');
    axis([0 1 0 1]);
    xlabel('Location x');
    ylabel('Concentration c');
    title(sprintf('Time step %d', i));
    pause(0.0005);
end

% ==========================================
% Local Functions for Thomas Algorithm
% ==========================================

function [P, G] = triDiagTrafo(R, A, B, C)
    % Initializes arrays for the new diagonal and right-hand side
    N_plus_1 = length(B);
    P = zeros(N_plus_1, 1); 
    G = zeros(N_plus_1 - 1, 1);
    
    % Initial row transformations
    P(1) = R(1) / B(1);
    G(1) = C(1) / B(1);
    
    % Transform the remaining rows to upper triangular form
    for i = 2:N_plus_1
        % A(1) mathematically points to a_2, requiring an i-1 shift
        denom = B(i) - A(i-1) * G(i-1); 
        P(i) = (R(i) - A(i-1) * P(i-1)) / denom;
        
        if i < N_plus_1
            G(i) = C(i) / denom;
        end
    end
end

function c_new = triDiagSolve(~, P, G)
    % Back-substitutes the upper triangular matrix to solve for c
    N_plus_1 = length(P);
    c_new = zeros(N_plus_1, 1);
    
    % Compute the right-most boundary node
    c_new(N_plus_1) = P(N_plus_1);
    
    % Compute the remaining interior nodes in reverse order
    for i = (N_plus_1 - 1):-1:1
        c_new(i) = P(i) - G(i) * c_new(i+1);
    end
end