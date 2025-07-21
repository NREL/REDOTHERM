function eta = OptimizeRedoxWrapper(X,S)
% This function is a wrapper for the optimization capability of REDOTHERM
% that runs the "Analyze_Redox_Cycle" function and only passes as an output
% the cycle efficiency (since this is the optimization variable)
Sol = Analyze_Redox_Cycle(X,S);
eta = -Sol.eta;
end