%
% This script can be used to plot the actual performance of the model,
% comparing its output to the real output one is trying to fit..
%
%
figure(1); clf;
modeloutput = model(input); % The model prediction is computed here..
P = plot(time, modeloutput,'r', time, real_output, 'k');
L = legend([P(1) P(2)], 'model', 'data');
