%
% This is the main fit routine..
%

global time input A EE H E;
global pars E_best p_best Nz Np;

% load_data;            <-- somewhere data should be made available

% Init values..
ZEROS = [0.1 20];
POLES = [10 250 250 300];
GAIN  = 3.4245;
DELAY = 0.001;
Nz    = length(ZEROS);
Np    = length(POLES);

%options = optimset('Display', 'on', 'TolX', 1e-9, 'TolFun', 1e-3, 'MaxFunEvals',10000, 'MaxIter', 10000);                       
%pars    = fminsearch(@cost, pars, options);

minZ = zeros(size(ZEROS));
minP = zeros(size(POLES));
maxZ = 10000*ones(size(ZEROS));
maxP = 20000*ones(size(POLES));
dpZ  = 20*ones(size(ZEROS));
dpP  = 20*ones(size(POLES));

pars  = [0.9366  172.2660    6.1414   -2.7866 GAIN DELAY ZEROS  POLES];
dpmax = [0.1    1     0.01      0.1 0.4     0.001   dpZ dpP];
p_min = [-5     0.    0.0001   -200. 0.0001 0    minZ minP];
p_max = [5      3000  10.       200. 1e3    0.002  maxZ maxP];

[E_best, p_best] = anneal(@cost, pars, dpmax, p_min, p_max, 100);

plot_model
