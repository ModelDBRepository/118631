function y = model(x)
%
% This assumes that the following data structures are declared as global vars: 
%
% time  --> a vector containing the time axis (e.g. time = 0:dt:Tmax, in [s]
% pars  --> a vector or parameters
% Nz    --> the number of parameters called zeros
% Np    --> the number of parameters called poles
%
% The function expects as input 
% a vector, same size of time, containing the input waveform to the 'model'
%

global time;                    % time axis is a global var..
global pars Nz Np;              % such as the model-fit parameters, and their number

N      = length(time);          % Total number of samples is determined..
dt     = time(2) - time(1);     % the sampling interval is determined.. [s]
omega  = 1./(N*dt);             % The sampling intervalin the frequency domain is determined [Hz]..
f      = omega*(0:(N-1))';      % The frequency-axis is generated here..

xx  = (x - mean(x))/1000.;      % Preprocessing of the input waveform
                                % i.e. removing its offset and rescaling
                                % This is of course equivalent to rescaling
                                % the zero-frequency gain of the 'model'
                                % It was used here for convenience.
                                
%---------------  STATIC NON-LINEARITY is defined here --------------------
xx  = (xx - pars(1));
xx  = pars(2) * (1./ (1 + exp(-xx * pars(3))));

%--------------- FILTER DESIGNS, WITH Nz ZEROS AND Np POLES ALL REAL ------
%
% Our convention here is that 'pars' contains a total of (Nz+Np+2) + 1 pars 
%
% The first 4 elements of 'pars' defines the sigmoid (see above)
% then the 5th element is the DC-gain, the 6th is the fixed time-delay
% From the 7th element on, Nz zeros are indicated and Np poles follow..
%

pars(6) = abs(mod(pars(6),0.1));    % Let's check that the delay is > 0
                                    % and that its step increase is 0.1 

NUM = 1;                            % I will use NUM and DEN to 'accumulate'
DEN = 1;                            % the filter, in the frequency domain
G   = 1;                            % See the definition in the Fourier..
FF  = zeros(size(f)) + sqrt(-1) * f;% FF = j * f

if (Nz > 0), for h=7:(Nz+7-1),          
        NUM = NUM .* (FF + pars(h));
        G   = G    * (+pars(h));    
end; end;
% i.e.
% NUM = (j*f + zero_1) * (j*f + zero_2) * ....
%
if (Np > 0), for h=(Nz+7-1+1):(Nz+Np+7-2+1),
        DEN = DEN .* (FF + pars(h));
        G   = G    * (+1./pars(h)); 
end; end;
% i.e.
% DEN = (j*f + pole_1) * (j*f + pole_2) * ....
%

PHASELAG = - f * pars(6);   % Constant time-delay <=> linear frequency phase
%
% i.e.     (j*f + zero_1) * (j*f + zero_2) * .... * (1/zero_1)*(1/zero_2)*..
%          ----------------------------------------------------------------
%          (j*f + pole_1) * (j*f + pole_2) * .... * (1/pole_1)*(1/pole_2)*..
%
T        = (pars(5)/G) * NUM ./ DEN;     
mT       = abs(T);       
pT       = angle(T);;

% In order to anti-transform and go back to the time-domain, I need to
% respect the Hilbert-symmetry for real signals..
k1 = 2;
k2 = N/2;
k3 = N;
mT((k1+k2):N) =  mT(k2:-1:k1);
pT((k1+k2):N) = -pT(k2:-1:k1);
PHASELAG((k1+k2):N) = - PHASELAG(k2:-1:k1);

% Finally the filter (in the Fourier Domain)
T = mT .* exp(sqrt(-1) * (pT+PHASELAG));
%%--------------------------------------------------------------------------------------------------------------

tmp = fft(xx) .* T;             % I/O Filtering in the frequency-domain is here!
y   = real(ifft(tmp));          % Let's go back to the time-domain

y   = y + pars(4);              % Last addition to the model (the DC level)
end