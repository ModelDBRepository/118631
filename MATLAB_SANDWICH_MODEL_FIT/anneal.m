function [E_best, p_best] = anneal(funfcn,p,dp_max,p_min,p_max,MAXITER)
%
%  ANNEAL Multidimensional constrained nonlinear minimization (Metropolis).
% [E_best p_best] = anneal(funfcn,p,dp_max,p_min,p_max,MAXITER)
%
% p      : Vector containing the current optimization parameters..
% dp_max : Maximal allowed modification for each parameter..
% p_max  : Allowed range for each parameter..
% p_min  : Allowed range for each parameter..
%
%
% p_best : Vector containing the 'best-ever' pars (i.e. corresponding to E_best)..
% E_best : Best ever value of the cost function
%
%   Reference: ???
%
%   Copyright 2006 Michele Giugliano, PhD
%   $Revision: 1.00 $  $Date: 2006/03/11 20:22:25 $

if nargin < 5, error('ANNEAL requires 5 input arguments'); end

global E_best p_best;

p_best = p;                % Vector containing the 'best-ever' pars (i.e. corresponding to E_best)..
t = cputime;
E      = feval(funfcn,p);  % Current value of the cost function (i.e. the current energy)..
t = cputime - t;
E_best = E;                % Best ever value of the cost function ('fake' initializatio to 1.e20)..
E_stuck_test    = 0.;      % Variable storing the value of the cost, when stuck in a local minimum..

disp(sprintf('\n\n\nE starting from %4.9f (%.0f ms to evaluate your function)\n', E,1000*t)); 

%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------
ii              = 1;       % Counter on the parameter number..
iter            = 2;       % Index counting iterations (it starts from 2)..
Npar            = prod(size(p));% Total no of *free* optimization parameters..
freeze          = 0;       % Boolean variable associated to high temperature "freezing".
Temp_init       = 4.;      % Initial value of the absolute 'temperature', for the annealing schedule.. 
Temp            = Temp_init;% Current value of the absolute 'temperature', for the annealing schedule..    
Temp_step       = 0.99;    % 'Temperature' decrease step factor, for the annealing schedule..
atten_init      = 1.;      % Initial value of the 'attenuation'..
atten           = atten_init;% Current value of the 'attenuation'..
atten_step      = 1.001;   % 'Attenuation' increase step factor, for the annealing schedule..

stuck_time      = 0;       % Variable counting the number of iterations, while stuck in a local minumum,.
max_stuck_times = 2000;    % Max number of times to be stucked in a local minimum..
restart         = 0;       % Variable storing the number of attempts at 'recoverying' from a local minimum..
max_restarts    = 10;      % Max number of times the process attempt at 'recoverying' from a local minimum..
max_iter        = 100000;  % Max overall number of iterations to be performed [UNUSED].
quake           = 100000;  % No of iterations after which an "earthquake" occurs.
quakes          = 2;       % Var storing the no of iteration before a minor 'earthquake'.

%MAXITER         = 100;     % Maximal number of allowed function iterations..
%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------

while(iter < MAXITER)                               % This is the main loop for the simulated annealing.
 iter = iter + 1;
 if (rem(iter, 100)==0), disp(sprintf('%d iterations so far..',iter)); end
%  disp(sprintf('%d iterations so far..',iter)); 
%---------
 if (~freeze)                                       % The 'temperature' is constantly decreased (for the next step)
  [Temp, atten] = annealing_schedule(Temp, Temp_step, atten, atten_step);
 end
%--------- 
 if (abs(E_stuck_test - E) < (1./1000.))           % If the parameters are stuck in a local minimum..
  stuck_time = stuck_time + 1;                      % then the change in the energy level is very small.
 else                                               % It is therefore useful to capture such a conditions,
  E_stuck_test = E;                                 % and when the number of stucked iterations is larger enough..
 end
%--------- 
 if (stuck_time == max_stuck_times)                 % It provides a 'temper. remescio' (jargon) to allow the minimum,
      Temp = Temp_init  * 0.5;                      % to be somehow overcome. Of course, at such a stage, all the
      atten= atten_init * 0.5;                      % (the 'atten' factor is also released a little bit).
      stuck_time = 0;      restart = restart + 1;   % counters must be reset to zero to start again the check.
 end
%---------
 if (restart == max_restarts)                       % However, after some time and after a given number of restarts
      Temp = Temp_init;       restart = 0;          % like those reported above, it is useful to have a large tempe-
      atten= atten_init;                            % rature outbreak as it was in the very beginning.
 end    
%---------    
 if (rem(quakes,quake)==0) 
  ii = fix(rand*(Npar-1)+1);                        % A parameter is chosen by chance..

  dp      =  dp_max(ii)* ( 1. - 2 * rand ) / atten;   % A random step (abs<dp_max(ii)) is attempted..      
  while (((p(ii)+dp) < p_min(ii)) | ((p(ii)+dp) > p_max(ii))) % Is the random step allowed ?
   dp      = dp_max(ii)* ( 1. - 2 * rand ) / atten;           % Otherwise, let's find another step.
  end
  p(ii)   = p(ii) + dp;                                       % When the step is ok, it is performed, so that  
  quakes = 0; 
 end    
%---------     
% Let's finally perform an evolution step..
 [E, E_best, p, p_best] = evolve(ii, funfcn, E, E_best, p, p_best, dp_max, p_min, p_max, Temp, atten);
 ii = max([rem(ii+1, Npar+1), 1]);
 quakes  = quakes + 1;
    
end
disp('Press enter to run it again.. CTRL-C to break!'); pause; [E_best, p_best] = anneal(funfcn,p_best,dp_max,p_min,p_max,MAXITER);
%------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------
function [Temp, atten] = annealing_schedule(Temp, Temp_step, atten, atten_step)   
  Temp = Temp * Temp_step;      % At each time step (i.e. each call) the temperature is lowered..
  atten= atten* atten_step;     % and the attenuation is increased.
%------------------------------------------------------------------------------------------------
   
%------------------------------------------------------------------------------------------------
function [A, B, C, D] = evolve(ii, funfcn, E, E_best, p, p_best, dp_max, p_min, p_max, Temp, atten)
%  if ((!divergence) && (ii==5))   % If no 'divergence' has been required, Tarp1 ..
%    ii = (ii+1) % Npar;           % The next parameter is welcome to change anyway..
A = E;     B = E_best;      C = p;      D = p_best;

 i = 1;
 dp      =  dp_max(ii)* ( 1. - 2 * rand ) / atten;   % A random step (abs<dp_max(ii)) is attempted..      
 while ((i<1000) | (((p(ii)+dp) < p_min(ii)) | ((p(ii)+dp) > p_max(ii)))) % Is the random step allowed ?
  dp      = dp_max(ii)* ( 1. - 2 * rand ) / atten;           % Otherwise, let's find another step.
  i = i + 1; 
 end
 p(ii)   = p(ii) + dp;                                       % When the step is ok, it is performed, so that  
 tmp     = feval(funfcn,p);
 if (isnan(tmp) | isinf(tmp)), return; end;
 
 dE      = E - tmp;                                          % the *descrease* in the cost is evaluated..
 if (dE > 0.)                                % If such a change in p(i) is convenient,
  E = E - dE;                                % we accept the new value for the p(i) and exit.
  if (E < E_best)                            % If this is even the best-ever score up to now,
   p_best = p;                               % we write down the corresponding parameters..
   disp(sprintf('****** E_best improved [%4.9f -> %4.9f]\n', E_best, E)); % and we let the user know it.
   E_best = E;                               % The new 'best-ever' is set here.
  end
 else                                        % When the step is not favorable (i.e. dE < 0)..
  if ( rand <= exp(dE/Temp) )                % let's *maybe* keep the change anyway (tossing a coin)..
   E = E - dE;                               % (this is the essence of the METROPOLIS algorithm)
  else
   p(ii) = p(ii) - dp;                       % Otherwise, discard the step.
  end
 end
A = E;     B = E_best;      C = p;      D = p_best;
%------------------------------------------------------------------------------------------------
