function out = cost(p)

global time input A pars;
pars = p;

N   = length(time);
out = (1./N) * sum( (( A - model(input) )).^2);
out = out + (1./N) * sum( (( diff(A,1) - diff(model(input),1) )).^2);
end