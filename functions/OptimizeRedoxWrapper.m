function eta = OptimizeRedoxWrapper(X,S)
% 
Sol = Analyze_Redox_Cycle(X,S);
eta = -Sol.eta;
end