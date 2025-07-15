function [X_sol,fvalue,exitflag] = OptimizeRedoxCycle(X0,S,opt_flag)
% Prepared by: Alon Lidor, NREL (alon.lidor@nrel.gov)
% This function is used to optimize the performance of a chemical looping
% redox cycle for H2O/CO2 splitting using redox-active material
% Input:
% X0 - an input vector consisting of: T_red, T_ox, p_red, p_ox, omega_red,
% omega_ox
% S - structure with all necessary parameters
% opt_flag - flag for choosing optimizer
lb = [S.T_red_min ; S.T_ox_min ; 0 ; 0 ; S.omega_red_min ; S.omega_ox_min];
ub = [S.T_red_max ; S.T_ox_max ; Inf ; Inf ; S.omega_red_max ; S.omega_ox_max];
Aeq = [0, 0, 1, 0, 0, 0 ; 
       0, 0, 0, 1, 0, 0];
beq = [S.p_red ; S.p_ox];
A = [-1, 1, 0, 0, 0, 0];
b = [0];
nvars = length(X0);
objective_function = @(X)OptimizeRedoxWrapper(X,S);

switch opt_flag
    case 1
        options = optimoptions('fmincon','PlotFcn',{'optimplotfunccount',...
            'optimplotfval','optimplotfirstorderopt'},'Algorithm','interior-point','HessianApproximation','lbfgs');
        [X_sol,fvalue,exitflag,output] = fmincon(objective_function,X0,A,b,Aeq,beq,lb,ub,[],options);
    case 2
        lb = [S.T_red_min ; S.T_ox_min ; S.p_red ; S.p_ox ; S.omega_red_min ; S.omega_ox_min];
        ub = [S.T_red_max ; S.T_ox_max ; S.p_red ; S.p_ox ; S.omega_red_max ; S.omega_ox_max];
        options = optimoptions('surrogateopt','Display','iter','ConstraintTolerance',1e-3,...
            'PlotFcn',{@optimplotfval,@optimplotx,@surrogateoptplot});
        secondary_options = optimoptions('fmincon','PlotFcn',{'optimplotfunccount',...
            'optimplotfval','optimplotfirstorderopt'},'Algorithm','interior-point','HessianApproximation','lbfgs');
        [X_sol_init,fvalue_init,exitflag_init,output_init,trials_init] = surrogateopt(objective_function,lb,ub,[],A,b,Aeq,beq,options);
        [X_sol,fvalue,exitflag,output] = fmincon(objective_function,X_sol_init,A,b,Aeq,beq,lb,ub,[],secondary_options);
    case 3
        options = optimoptions('patternsearch','Algorithm','nups','Display','iter','Cache','on',...
            'MaxFunctionEvaluations',5000,'UseParallel',false,'UseCompleteSearch',true,'ScaleMesh',false,...
            'UseCompletePoll',true,'PlotFcn',{'psplotbestf','psplotbestx'});
        [X_sol,fvalue,exitflag,output] = patternsearch(objective_function,X0,A,b,Aeq,beq,lb,ub,[],options);
    case 4
        hybrdidopts = optimoptions('fmincon','OptimalityTolerance',1e-9,...
            'PlotFcn',{'optimplotfval','optimplotx'});
        options = optimoptions('ga','HybridFcn',{@fmincon,hybrdidopts},'Display','iter','PlotFcn',{@gaplotbestf},...
            'ConstraintTolerance',1e-6,'FunctionTolerance',1e-6,'UseVectorized',false);
        % options = optimoptions('ga','Display','iter','PlotFcn',{'gaplotdistance','gaplotbestf'},...
        %     'ConstraintTolerance',1e-6,'FunctionTolerance',1e-6,'UseVectorized',false);
        [X_sol,fvalue,exitflag,output,population,scores] = ga(objective_function,nvars,A,b,Aeq,beq,lb,ub,[],options);
    case 5
        options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
            'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt},'ConstraintTolerance',1e-6,...
            'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
        % Problem structure
        problem = createOptimProblem('fmincon','x0',X0,'objective',objective_function,...
            'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
        % Solve
        gs = GlobalSearch;
        gs.Display = 'iter';
%         gs.OutputFcn = @outputfun;
        gs.PlotFcn = {@gsplotbestf,@gsplotfunccount};
        gs.NumTrialPoints = 2000;
        gs.NumStageOnePoints = 250;
        % gs.MaxTime = 300; % s (5min)
        gs.XTolerance = 1e-6;
        gs.FunctionTolerance = 1e-6;
        gs.StartPointsToRun = 'bounds-ineqs';
%         gs;
        % Run GlobalSearch
        [X_sol,favlue,exitflag,output,solutions] = run(gs,problem);
    case 6
        % Local solver options
        options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
            'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt},'ConstraintTolerance',1e-6,...
            'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
        % Problem structure
        problem = createOptimProblem('fmincon','x0',X0,'objective',objective_function,...
            'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
        % Solve
        gs = GlobalSearch;
        gs.Display = 'iter';
%         gs.OutputFcn = @outputfun;
%         gs.PlotFcn = {@gsplotbestf,@gsplotfunccount};
        gs.NumTrialPoints = 1e5;
        % gs.MaxTime = 300; % s (5min)
        gs.XTolerance = 1e-6;
        gs.FunctionTolerance = 1e-6;
        gs.StartPointsToRun = 'bounds-ineqs';
%             gs;
        % MultiStart
        ms = MultiStart(gs);
        ms.UseParallel = true; % Attempet to use parallel computing
        % ms.PlotFcn = {@gsplotbestf,@gsplotfunccount};
        p = gcp('nocreate');
        if isempty(p)
            % There is no parallel pool
            parpool;
        end
%         ms;
        stpoints = RandomStartPointSet('NumStartPoints',10);    % Random starting points
        % Run the MultiStart solver
        disp('Optimizing...');
        tic;
        [X_sol,fvalue,exitflag,output,solutions] = run(ms,problem,stpoints);
        delete(gcp);
    case 7
        p = gcp('nocreate');
        if isempty(p)
            % There is no parallel pool
            parpool;
        end
        options = optimoptions('patternsearch','Algorithm','nups','Display','iter','Cache','on',...
            'MaxFunctionEvaluations',5000,'UseParallel',true,'UseCompleteSearch',true,...
            'UseCompletePoll',true,'PlotFcn',{'psplotbestf','psplotbestx'});
        [X_sol,fvalue,exitflag,output] = patternsearch(objective_function,X0,A,b,Aeq,beq,lb,ub,[],options);
        delete(gcp);
    case 8
        % Local solver options
        options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
            'PlotFcn',{@optimplotfval,@optimplotx,@optimplotfirstorderopt},'ConstraintTolerance',1e-6,...
            'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
        % Problem structure
        problem = createOptimProblem('fmincon','x0',X0,'objective',objective_function,...
            'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
        % Solve
        gs = GlobalSearch;
        gs.Display = 'iter';
%         gs.OutputFcn = @outputfun;
%         gs.PlotFcn = {@gsplotbestf,@gsplotfunccount};
        gs.NumTrialPoints = 1e5;
        % gs.MaxTime = 300; % s (5min)
        gs.XTolerance = 1e-6;
        gs.FunctionTolerance = 1e-6;
        gs.StartPointsToRun = 'bounds-ineqs';
%             gs;
        % MultiStart
        ms = MultiStart(gs);
        % ms.PlotFcn = {@gsplotbestf,@gsplotfunccount};
%         ms;
        stpoints = RandomStartPointSet('NumStartPoints',10);    % Random starting points
        % Run the MultiStart solver
        disp('Optimizing...');
        tic;
        [X_sol,fvalue,exitflag,output,solutions] = run(ms,problem,stpoints);
end

end