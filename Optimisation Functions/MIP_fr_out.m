% MIP_fr_out
% ==========

%  Purpose
%  -------
%  Maximises benefits suject to a budget constraint, by offering flat rate 
%  prices to farmers for each outcome arising from farmer choice over a set 
%  of possible farm land use change options. Outcomes can be for
%  environmental quantities (fr_env) or for ecosystem service values
%  (fr_es).

%  Inputs
%  ------
%    b                  [num_farmers x num_options] benefits from each option 
%    c                  [num_farmers x num_options] costs of each option 
%    q                  [num_farmers x num_env_out x num_options] quantities of env outs from each option  
%    budget             Maximum spend   
%    lu_data            [num_farmers x num_options] indicator of which options refer to the land use in this observation (all for cell, ag/grass for parcel) 
%    warm_start_prices  [1 x num_env_out] vector of best starting prices
%    warm_start_uptake  [1 x num_farmers*num_options] vector of uptake at warm_start_prices 
%    prices_lb          [1 x num_env_out] vector of lower bount on prices
%    prices_ub          [1 x num_env_out] vector of upper bount on prices
%    cplex_options      structure with cplex options
%                         o time: secs to allow cplex to run 
%                         o logs: folder in which to find warmstarts and
%                                 write node logs

function [prices, uptake, fval, exitflag, exitmsg] =  MIP_fr_out(b, c, q, budget, lu_data, warm_start_prices, warm_start_uptake, prices_lb, prices_ub, cnst_data, cnst_target, byparcel, cplex_options)
    
    % 1. Initialise
    % =============
    
    % 1.1 Cplex Object
    % ----------------
    if ismac
        addpath(genpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx'))
    elseif ispc
        addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64'))
    end
    
    cplex = Cplex('elms_lp');
    cplex.Model.sense = 'maximize';
    cplex.Param.emphasis.mip.Cur = 0; % balanced
	% cplex.Param.emphasis.mip.Cur = 1; % emphasise feasibility
    cplex.Param.mip.strategy.search.Cur = 2;
    cplex.Param.parallel.Cur = 1;
    cplex.Param.mip.tolerances.integrality.Cur = 0;
    cplex.Param.mip.tolerances.mipgap.Cur = 0;
    cplex.Param.timelimit.Cur = cplex_options.time;
    cplex.Param.workdir.Cur = cplex_options.logs;
    %  cplex.DisplayFunc = '';    
    
    % 1.2 Constants
    % -------------
    num_farmers = size(q, 1);
    num_env_out = size(q, 2);
    num_options = size(q, 3);
    
    % 1.3 Reduce to Cells where could have positive surplus given prices
    % ------------------------------------------------------------------
    % max_surplus = fcn_find_max_surplus(q, c, prices_lb, prices_ub);
    surplus = fcn_find_max_surplus(q, c, prices_lb, prices_ub);
    max_surplus = max(surplus, [], 2);
    
    drop_cell_ind = (max_surplus < 0);       
    if any(drop_cell_ind)
        %  Remove cells where could no possible price would enduce participation
        b(drop_cell_ind, :)    = [];
        c(drop_cell_ind, :)    = [];
        q(drop_cell_ind, :, :) = [];
        lu_data(drop_cell_ind, :)    = [];
        cnst_data(:,:,drop_cell_ind) = [];
        surplus(drop_cell_ind, :)  = [];
        max_surplus(drop_cell_ind) = [];
        warm_start_uptake(repelem(drop_cell_ind, num_options, 1)') = [];
        num_farmers = size(q, 1);
    end
    R = repelem(budget*1.01, num_options * num_farmers, 1);
    S = budget*1.01;
    
    % R = repelem(max_surplus+1, num_options, 1);
    % S = max_surplus+1;
    
    % Store transposed matrices
    bt = b';
    ct = c';    
    
   
    % 2. Choice Variables
    % ===================
    %
    %    Prices     Farm Utility       Activity Choice
    %   --------   --------------   ---------------------  
    %                                Farm1   ...  FarmN 
    % [1 x Nprice]    [1 x N]      [1 x Nopt]   [1 x Nopt]
    %      >0            >0           0-1          0-1 
    
    
    % 3. Bounds for Choice Variables
    % ==============================
    farm_option_ub = lu_data';
    farm_option_ub = farm_option_ub(:);
    lb    = [prices_lb; sparse(num_farmers + num_options * num_farmers, 1)];
    ub    = [prices_ub; repelem(budget, num_farmers, 1); farm_option_ub];     
   
    
    % 4. Cost vector
    % ==============    
    f_p = sparse(num_env_out, 1);
    f_u = sparse(num_farmers, 1);
    f_x = bt(:);
    f = [f_p; f_u; f_x];
    clear f_p f_u f_x

    
    % 5. Variable types
    % =================
    ctype = [repmat('C',1,(num_env_out + num_farmers)), repmat('B', 1, (num_options * num_farmers))];

    
    cplex.addCols(f, [], lb, ub, ctype);  
    
    clear f lb ub ctype
    
    
    % 6. Inequality Constraints
    % =========================

    % 6.1 Unit Demand
    % ---------------
    %  Farm can select at most one of the options.
    %
    %     sum_f(x_of) <= 1   (f = 1...N)      (Eq 12)
    %    
    A_p = sparse(num_farmers, num_env_out);
    A_u = sparse(num_farmers, num_farmers);
    A_x = repelem(speye(num_farmers), 1, num_options);
    A   = [A_p, A_u, A_x];
    Bl  = sparse(num_farmers,1);
    Bu  = ones(num_farmers,1);
    if ~isempty(A),cplex.addRows(Bl, A, Bu); end
    clear A A_p A_x Bl Bu

    % 6.2 Utility Maximisation: Lower Bound
    % -------------------------------------
    %  Farm utility must be >= utility of option offering the most profit
    %  at current prices
    %
    %    u_f >= sum_e(p_e.q_eof) - c_of) (f = 1...N; o = 1...No)   (Eq 13)
    %    
    A_p = permute(q, [3, 2, 1]);
    A_p = reshape(permute(A_p, [1, 3, 2]), [num_farmers*num_options num_env_out]);   
    A_u = -repelem(speye(num_farmers), num_options, 1);
    A_x = sparse(num_farmers*num_options, num_farmers*num_options);
    A   = [A_p, A_u, A_x];
    Bl  = ones(num_farmers*num_options, 1) * -Inf;
    Bu  = ct(:);
    if ~isempty(A),cplex.addRows(Bl, A, Bu); end
    clear A A_u A_x Bl Bu    
    
    
    % 6.3 Utility Maximisation: Upper Bound
    % -------------------------------------
    %  Farm utility <= utility of option chosen by farmer, tying the choice
    %  of option to utility maximising choice from (Eq 12).
    %
    %    u_f <= sum_e(p_e.q_eof) - sum_oo(c_oof.x_oof) + (1 - x_of)R  (f = 1...N; o = 1...No)   (Eq 14)
    %   
    %  Here R is a constant that is sufficiently high to ensure that if
    %  option o is not chosen (x_of = 0) that the rhs holds trivially. Only
    %  when the option is chosen (x_of = 1) will this upper bound be placed
    %  on farmer utility. The only way the optimiser can get (6.2) and
    %  (6.3) to hold simultaneously is if the chosen option is also the one
    %  that gives the highest utility of all options. Alternatively, if all
    %  options give negative utility then u_f = 0 on account of positivity
    %  constraint, in which case all x_of must be zero and no option is
    %  chosen.    
    % R   = budget .* 1.01;
    A_p = -A_p;   
    A_u = repelem(speye(num_farmers), num_options, 1);
    A_x = num2cell(c,2);
    A_x = repelem(sparse(blkdiag(A_x{:})), num_options, 1);
    A_x = A_x + R.*speye(num_farmers*num_options);
    A   = [A_p, A_u, A_x];
    B   = ones(num_farmers*num_options, 1);
    Bl  = B * -Inf;
    Bu  = B .* R;
    if ~isempty(A),cplex.addRows(Bl, A, Bu); end    
    clear A_p A_u A_x c_rf A B Bl Bu   

    
    % 6.4 Budget Constraint
    % ---------------------
    %  Amount paid to farmers cannot exceed budget.
    %
    %    sum_f(u_f) + sum_f(sum_o(c_of.x_of)) <= M   (Eq 17)
    %   
    A_p = sparse(1, num_env_out);
    A_u = ones(1, num_farmers);
    A_x = ct(:)';
    A   = [A_p, A_u, A_x];
    Bl  = 0;
    Bu  = budget;
    if ~isempty(A),cplex.addRows(Bl, A, Bu); end        
    clear A_p A_u A_x A Bl Bu  
    
    % 6.5 Biodiversity Constraints
    % ----------------------------
    % Increases in cell counts of species within each functional group must
    % by greater than the target
    if any(cnst_target)
        num_groups = length(cnst_target); 
        cnst_data = reshape(cnst_data, num_groups, []); % Reshape so have blocks of option to species group outcomes for each cell
        A_p = sparse(num_groups, num_env_out);
        A_u = sparse(num_groups, num_farmers);
        A_x = cnst_data;
        A  = [A_p, A_u, A_x];
        Bl = cnst_target;
        Bu = inf(num_groups,1);
        if ~isempty(A),cplex.addRows(Bl, A, Bu); end
        clear A A_p A_x Bl Bu    
    end       
    
    % 6.6 Additional Cuts
    % -------------------
    %  Surplus must be zero if no option chosen
    %
    %    u_f <= S * sum_o(x_of)   (f = 1...N)   (Eq 16)
    %   
    % S   = budget * 1.01;
    A_p = sparse(num_farmers, num_env_out);
    A_u = speye(num_farmers);
    A_x = -repelem(speye(num_farmers), 1, num_options) .* S;
    A   = [A_p, A_u, A_x];
    Bl  = ones(num_farmers, 1) * -Inf;
    Bu  = zeros(num_farmers, 1);
    if ~isempty(A),cplex.addRows(Bl, A, Bu); end                
    clear A_p A_u A_x A Bl Bu

%     %  Surplus must be > 1 if option chosen
%     %
%     %    u_f >= sum_o(x_of)   (f = 1...N)   (Eq 18)
%     %   
%     A_p = sparse(num_farmers, num_env_out);
%     A_u = -speye(num_farmers);
%     A_x = -repelem(speye(num_farmers), 1, num_options);
%     A   = [A_p, A_u, A_x];
%     Bl  = ones(num_farmers, 1) * -Inf;
%     Bu  = zeros(num_farmers, 1);
%     if ~isempty(A),cplex.addRows(Bl, A, Bu); end    
%     clear A_p A_u A_x A Bl Bu  




    % 7. CPLEX call
    % =============
    
    % 7.1 Warm start
    % --------------
    if ~isempty(warm_start_prices)
        sln = [warm_start_prices; warm_start_uptake]';
        idx = [1:num_env_out num_env_out+num_farmers+1:num_env_out+num_farmers+num_options*num_farmers];
        idx = idx - 1;            
        filename = [cplex_options.logs 'warmstart.mst'];
        probname = 'elms_lp';
        fcn_write_warmstart(sln', idx', filename, probname);
        cplex.readMipStart(filename);
        cplex.Param.mip.limits.repairtries.Cur = 500; % set the number of MIPstart repair attempts to 10
    end
    
    % 7.2 Solve
    % ---------
    tic    
    cplex.solve();
    toc   
    
    % 7.3 Collect Outputs
    % -------------------
    x        = cplex.Solution.x;
    prices   = x(1:num_env_out)';
    uptake   = round(reshape(x((num_env_out + num_farmers + 1):end), num_options, num_farmers)'); 
    fval     = cplex.Solution.objval; 
    exitflag = cplex.Solution.status;
    exitmsg  = cplex.Solution.statusstring;
    
end
