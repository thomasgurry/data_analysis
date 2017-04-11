%% Input template for a reaction system for use with the functions 
%% LNA and MRE.  The system is a single gene expression system in
%% this case.

% DNA ---[kd]---> DNA + r
%  r  --[gr*r]-->    0
%  r  --[kp*r]-->  r + p
%  p  --[gp*p]-->    0

global S macroF

% Define parameters

k1 = 1;
k2 = 0.1;
k3 = 1;
k4 = 0.01;

param = {'k1','k2','k3','k4'};

% Define the stoichiometry matrix

S = [1 -1 0 0 ; 0 0 1 -1];

% Define macroscopic transition rates (macroF).  Here kd = 1, gr = 0.01, 
% kp = 1 and gp = 0.01.  The variable names of the species must be in the
% form 'x1x' for species 1 (here x1x = r), 'x2x' for species 2 (here x2x = p), etc.

Sdim = size(S);
nreactions = Sdim(2);
macroF = sym(zeros(1,nreactions));
macroF(1) = 'k1';
macroF(2) = 'k2*x1x';
macroF(3) = 'k3*x1x';
macroF(4) = 'k4*x2x';

% Define the initial conditions

xinit = [0 0];