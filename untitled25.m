% test_3chamber_admittance.m
% Isolated verification of 3-chamber cochlear admittance mapping
clc; clear;

% 1. Arbitrary partition masses
m1 = 1.0; % Basilar Membrane mass
m2 = 0.5; % Tectorial Membrane mass
mu = m1 / m2;

% 2. The Native Difference Matrix D (From fdm24.m)
% Maps fluid pressure to mechanical degrees of freedom (P_mech = D * P_fluid)
D = [ 1,  0, -1;  % BM feels P_SV - P_ST
      0, -1,  1]; % TM feels P_ST - P_SS

% =========================================================================
% TEST 1: The tdm26.m Baseline (The Hydraulic Lock)
% =========================================================================
% In tdm26.m, fold_p defines Hair Bundle acceleration as:
% a(i2) = s1/m1 + s2/m2 
% This means the inverse mass matrix (M_inv) is NOT diagonal.
M_inv_tdm26 = [1/m1,    0;
               1/m1, 1/m2];

% Calculate resulting fluid admittance: Y = D^T * M_inv * D
Y_tdm26 = D' * M_inv_tdm26 * D;

disp('--- TEST 1: Current tdm26.m Admittance ---');
disp('Y_fluid = '); disp(Y_tdm26);
disp('Is the matrix symmetric? (Maxwell-Betti Reciprocity)');
if issymmetric(Y_tdm26), disp('TRUE'); else, disp('FALSE - HYDRAULIC LOCK'); end
disp(' ');

% =========================================================================
% TEST 2: The Strictly Diagonal Solution
% =========================================================================
% If fold_p is updated so kinematics are strictly independent:
% a(i2) = s2/m2
M_inv_sym = [1/m1,    0;
                0, 1/m2];

Y_sym = D' * M_inv_sym * D;

disp('--- TEST 2: Corrected Symmetric Admittance ---');
disp('Y_fluid = '); disp(Y_sym);
disp('Is the matrix symmetric?');
if issymmetric(Y_sym), disp('TRUE - STABLE'); else, disp('FALSE'); end

% Normalize by m1 to match the time-domain spatial diagonals (abmom)
Y_norm = m1 * Y_sym;
disp(' ');
disp('--- Normalized Y_norm (The 3x3 block for cochlea.m) ---');
disp(Y_norm);
disp('Do the columns and rows sum to 0? (Volume Conservation)');
disp(sum(Y_norm)); 

% =========================================================================
% TEST 3: Forcing Vector Routing
% =========================================================================
% Fluid forces are routed by Q_fluid = D^T * M_inv * S_mech
% S_norm is Q_fluid normalized by m1.
disp(' ');
disp('--- Normalized Forcing Vector S_norm (for xpnd_q) ---');
disp('Q_SV =  s1');
disp('Q_SS = -mu * s2');
disp('Q_ST = -s1 + mu * s2');