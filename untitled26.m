% test_fluid_sandbox.m
% Isolated test of the fully coupled 3-chamber fluid matrix
clc; clear; close all;

N = 100; % Spatial nodes
M = 3;   % Chambers
disp('Building 3-chamber fully coupled fluid matrix...');

% 1. Dummy physical parameters (constant for simplicity)
aflom = 0.1 * ones(N,1); % Acoustic mass
abmom = 1.0 * ones(N,1); % BM mass
mu    = 0.5 * ones(N,1); % Mass ratio (m1/m2)
alfx  = 0.5;             % ME coupling

a1 = zeros(N, M*M); a2 = zeros(N, M*M); a3 = zeros(N, M*M);
kk = 2:(N-1);

% 2. Build the Spatial Laplacians
for k = 1:M
    diag_idx = k + (k-1)*M;
    a1(kk, diag_idx) = -aflom(kk-1) ./ abmom(kk);
    a3(kk, diag_idx) = -aflom(kk)   ./ abmom(kk);
    a2(kk, diag_idx) = -a1(kk, diag_idx) - a3(kk, diag_idx); 
end

% 3. Superimpose the Symmetric Y_norm (The cross-chamber fluid coupling)
a2(kk, 1) = a2(kk, 1) + 1;
a2(kk, 3) = a2(kk, 3) - 1;
a2(kk, 5) = a2(kk, 5) + mu(kk);
a2(kk, 6) = a2(kk, 6) - mu(kk);
a2(kk, 7) = a2(kk, 7) - 1;
a2(kk, 8) = a2(kk, 8) - mu(kk);
a2(kk, 9) = a2(kk, 9) + 1 + mu(kk);

% 4. Basal BC (Node 1)
for k = 1:M
    diag_idx = k + (k-1)*M;
    a3(1, diag_idx) = -aflom(1) ./ abmom(1);
    a2(1, diag_idx) = -a3(1, diag_idx);
end
a2(1, 1) = a2(1, 1) + 1; a2(1, 3) = a2(1, 3) - 1;
a2(1, 5) = a2(1, 5) + mu(1); a2(1, 6) = a2(1, 6) - mu(1);
a2(1, 7) = a2(1, 7) - 1; a2(1, 8) = a2(1, 8) - mu(1);
a2(1, 9) = a2(1, 9) + 1 + mu(1);

a2(1, 1) = a2(1, 1) + alfx; % ME Drive
a3(1, 9) = 0; a2(1, 7:9) = 0; a2(1, 9) = 1; % Round Window Clamp P_ST=0

% 5. Apical BC (Node N) - Helicotrema
for k = 1:M
    diag_idx = k + (k-1)*M;
    a1(N, diag_idx) = -1;
    a2(N, diag_idx) =  1;
end

% Assemble Sparse Matrix
AA = sandbox_xpnd_a(a1, a2, a3, M, N);
disp('Matrix Condition Number (Lower is better, >1e15 is singular):');
disp(condest(AA));

% 6. Apply a static Stapes Drive
Q = zeros(N*M, 1);
Q(1) = 1.0; % +Drive to SV
Q(3) = -1.0; % -Return from ST
Q(3) = 0; % Force RW boundary RHS to 0 to match clamp

% Solve for Pressures
P = AA \ Q;
P_SV = P(1:3:end);
P_SS = P(2:3:end);
P_ST = P(3:3:end);

% Plot the hydrodynamics
figure;
plot(1:N, P_SV, 'r', 'LineWidth', 2); hold on;
plot(1:N, P_SS, 'g', 'LineWidth', 2);
plot(1:N, P_ST, 'b', 'LineWidth', 2);
legend('Scala Vestibuli', 'Spiral Sulcus', 'Scala Tympani');
title('Static Fluid Pressure Gradient (Sandbox)');
xlabel('Spatial Node'); ylabel('Pressure');
grid on;

% ---- Helper Function ----
function aa=sandbox_xpnd_a(a1,a2,a3,m,n)
    nm = n*m; nd = 1+2*m; ad = zeros(nm,nd); dd = zeros(nm,nd);
    for j=1:m
        jj = (1:m)+(j-1)*m; kk_vec = (j:m:nm)';
        AAA_j = [a1(:, jj), a2(:, jj), a3(:, jj)];
        src_indices = (1:nd) + (j-1); 
        ad(kk_vec, :) = AAA_j(:, src_indices);
    end
    dd(:, 1+m) = ad(:, 1+m);
    for j=1:m
        idx_low = 1+m-j; idx_high = 1+m+j;
        dd(1:(nm-j), idx_low) = ad((1+j):nm, idx_low);
        dd((1+j):nm, idx_high) = ad(1:(nm-j), idx_high);
    end
    aa = spdiags(dd,-m:m,nm,nm);
end