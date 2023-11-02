% MIP_fr_act
% ==========

%  Purpose
%  -------
%  Maximises benefits suject to a budget constraint, by offering flat rate 
%  prices to farmers for each different land use change option (fr_es).

%  Inputs
%  ------
%    b                  [num_farmers x num_options] benefits from each option 
%    c                  [num_farmers x num_options] costs of each option 
%    q                  [num_farmers x num_options x num_options] quantities of env outs from each option  
%    budget             Maximum spend   
%    lu_data            [num_farmers x num_options] indicator of which options refer to the land use in this observation (all for cell, ag/grass for parcel) 
%    warm_start_prices  [1 x num_env_out] vector of best starting prices
%    warm_start_uptake  [1 x num_farmers*num_options] vector of uptake at warm_start_prices 
%    prices_lb          [1 x num_options] vector of lower bound on prices
%    prices_ub          [1 x num_options] vector of upper bound on prices
%    cplex_options      structure with cplex options
%                         o time: secs to allow cplex to run 
%                         o logs: folder in which to find warmstarts and
%                                 write node logs

function [prices, uptake, fval, exitflag, exitmsg] = MIP_fr_act(b, c, q, budget, lu_data, warm_start_prices, warm_start_uptake, prices_lb, prices_ub, cnst_data, cnst_target, cplex_options)


    % 1. Initialise
    % =============
    
    % 1.1. Add the cplex path into matlab
    if ismac
        addpath(genpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx'))
    elseif ispc
        addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64'))
    end
      
    % 1.2. CLPEX optimisation object
    % ------------------------------
    cplex = Cplex('elms_lp');
    cplex.Model.sense = 'maximize';
    cplex.Param.emphasis.mip.Cur = 0; % balanced
	% cplex.Param.emphasis.mip.Cur = 1; % emphasise feasibility
    cplex.Param.mip.strategy.search.Cur = 2;
    cplex.Param.parallel.Cur = 1;
    cplex.Param.mip.tolerances.integrality.Cur = 0;
    cplex.Param.mip.tolerances.mipgap.Cur = 0;
    cplex.Param.mip.tolerances.absmipgap.Cur = 0;
    cplex.Param.timelimit.Cur = cplex_options.time;
    cplex.Param.workdir.Cur   = cplex_options.logs;
    %  cplex.DisplayFunc = '';   
    
    % 1.2 Constants
    % -------------
    num_farmers = size(b, 1);
    num_options = size(b, 2);
    num_x = num_options + num_farmers + num_options * num_farmers;
    
    % Store transposed matrices
    bt = b';
    ct = c';
    qt = q';

    % set the parameters for the big-M formulation
    qf = max(q, [], 2);
    C_ha = (c ./ q);
    Ro = sparse(1, num_options);
    for col = 1: num_options
        C_col = C_ha(:,col);
        Ro(col) = max(C_col(~isinf(C_col)));
    end
    R = max(Ro(~isinf(Ro))) * 10;
    Rf = R * qf;
    S = max(Rf);
    T = 1;

    
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
    lb = [prices_lb; sparse(num_farmers + num_options * num_farmers, 1)];
    ub = [prices_ub; repelem(budget, num_farmers, 1); farm_option_ub];     
    

    % 4. Cost vector
    % ==============    
    f_p = sparse(num_options, 1);
    f_u = sparse(num_farmers, 1);
    f_x = bt(:);
    f = [f_p; f_u; f_x];
    clear f_p f_u f_x
   
    % 5. Variable types
    % =================
    ctype = [repmat('C',1,(num_options + num_farmers)), repmat('B', 1, (num_options * num_farmers))];
    
    % Add to cplex class
    cplex.addCols(f, [], lb, ub, ctype);

    
    % 6. Inequality Constraints
    % =========================

    % 6.1 Unit Demand
    % ---------------
    Aineq1_p = sparse(num_farmers, num_options);
    Aineq1_u = sparse(num_farmers, num_farmers);
    Aineq1_d = kron(speye(num_farmers), ones(1, num_options));
    Aineq1 = [Aineq1_p, Aineq1_u, Aineq1_d];

    Bineq1 = ones(num_farmers, 1);
    B1_lb  = zeros(num_farmers, 1);
    B1_ub  = Bineq1;
    
    clear Aineq1_p Aineq1_u Aineq1_d

    % 2nd inequality
    Aineq2_p = [];
    for i = 1:num_farmers
        Aineq2_p = [Aineq2_p; spdiags(q(i, :)', 0, num_options, num_options)];
    end
    Aineq2_u = -repelem(speye(num_farmers), num_options, 1);
    Aineq2_d = sparse(num_options * num_farmers,num_options * num_farmers);
    Aineq2 = [Aineq2_p, Aineq2_u, Aineq2_d];

    Bineq2 = reshape(c', [num_options*num_farmers, 1]);
    B2_lb = ones(length(Bineq2), 1) * -Inf;
    B2_ub = Bineq2;

    clear Aineq2_p Aineq2_u Aineq2_d

    % 3rd inequality
    Aineq3_p = [];
    for i = 1:num_farmers
        Aineq3_p = [Aineq3_p; -spdiags(q(i,:)', 0, num_options, num_options)];
    end
    Aineq3_u = repelem(speye(num_farmers), num_options, 1);

    eye_temp = speye(num_farmers);
    Aineq3_d = [];
    for i = 1:num_farmers
        Aineq3_d = [Aineq3_d, kron(c(i, :), eye_temp(i, :)')];
    end
    Aineq3_d = repelem(Aineq3_d, num_options, 1);

    Rf_rep = repmat(Rf, 1, num_options)';
    Rf_rep = Rf_rep(:);
    Aineq3_d = Aineq3_d + (Rf_rep.*speye(num_farmers*num_options));
    Aineq3 = [Aineq3_p, Aineq3_u, Aineq3_d];

    Bineq3 = Rf_rep;
    B3_lb = ones(length(Bineq3), 1) * -Inf;
    B3_ub = Bineq3;

    clear Aineq3_p Aineq3_u Aineq3_d

    % 4th inequality
    Aineq4_p = -speye(num_options);
    Aineq4_u = sparse(num_options, num_farmers);
    Aineq4_d = sparse(num_options, num_farmers * num_options);
    Aineq4 = [Aineq4_p, Aineq4_u, Aineq4_d];

    Bineq4 = sparse(num_options, 1);
    B4_lb = ones(num_options, 1) * -Inf;
    B4_ub = Bineq4;

    clear Aineq4_p Aineq4_u Aineq4_d

    % 5th inequality
    Aineq5_p = sparse(num_farmers, num_options);
    Aineq5_u = -speye(num_farmers);
    Aineq5_d = sparse(num_farmers, num_farmers * num_options);
    Aineq5 = [Aineq5_p, Aineq5_u, Aineq5_d];

    Bineq5 = sparse(num_farmers, 1);
    B5_lb = ones(num_farmers, 1) * -Inf;
    B5_ub = Bineq5;

    clear Aineq5_p Aineq5_u Aineq5_d

    % 6th inequality constraint: budget
    Aineq6_p = sparse(1, num_options);
    Aineq6_u = ones(1, num_farmers);
    Aineq6_d = ct(:)';
    Aineq6 = [Aineq6_p, Aineq6_u, Aineq6_d];

    Bineq6 = budget;
    B6_lb = 0;
    B6_ub = budget;
%     B6_lb = budget;
%     B6_ub = +inf;
    
    clear Aineq6_p Aineq6_u Aineq6_d

   
    % 8th inequality constraint: x = 1 when Uf>0
    Aineq8_p = zeros(num_farmers, num_options);
    Aineq8_u = speye(num_farmers);
    Aineq8_d = -kron(speye(num_farmers), ones(1, num_options)) .* S;
    Aineq8 = [Aineq8_p, Aineq8_u, Aineq8_d];

    Bineq8 = zeros(num_farmers, 1);
    B8_lb = ones(num_farmers, 1) * -Inf;
    B8_ub = Bineq8;
    
    clear Aineq8_p Aineq8_u Aineq8_d

    % 9th constraint (Utility must be slightly greater than 0 for having an uptake)
    Aineq9_p = zeros(num_farmers, num_options);
    Aineq9_u = -speye(num_farmers) .* T;
    Aineq9_d = kron(speye(num_farmers), ones(1, num_options));
    Aineq9 = [Aineq9_p, Aineq9_u, Aineq9_d];

    Bineq9 = zeros(num_farmers, 1);
    B9_lb = ones(num_farmers, 1) * -Inf;
    B9_ub = Bineq9;
    
    clear Aineq9_p Aineq9_u Aineq9_d
    
    % 10th inequality constraint: p_o=0 when sum over f(x_of)=0
    Aineq10_p = speye(num_options);
    Aineq10_u = zeros(num_options, num_farmers);
    Aineq10_d = -kron(ones(1, num_farmers), eye(num_options)) .* S;
    Aineq10 = [Aineq10_p, Aineq10_u, Aineq10_d];

    Bineq10 = zeros(num_options, 1);
    B10_lb = ones(num_options, 1) * -Inf;
    B10_ub = Bineq10;
    
    clear Aineq10_p Aineq10_u Aineq10_d
    
    % Combine all inequalities into one matrix
    A = [Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6; Aineq8; Aineq9; Aineq10];
    B = [Bineq1; Bineq2; Bineq3; Bineq4; Bineq5; Bineq6; Bineq8; Bineq9; Bineq10];
    B_lb = [B1_lb; B2_lb; B3_lb; B4_lb; B5_lb; B6_lb; B8_lb; B9_lb; B10_lb];
    B_ub = [B1_ub; B2_ub; B3_ub; B4_ub; B5_ub; B6_ub; B8_ub; B9_ub; B10_ub];
    cplex.addRows(B_lb, A, B_ub);

    clear Aineq1 Aineq2 Aineq3 Aineq4 Aineq5 Aineq6 Aineq8 Aineq9 Aineq10
    clear Bineq1 Bineq2 Bineq3 Bineq4 Bineq5 Bineq6 Bineq8 Bineq9 Bineq10
    clear B1_lb B2_lb B3_lb B4_lb B5_lb B6_lb B8_lb B9_lb B10_lb
    clear B1_ub B2_ub B3_ub B4_ub B5_ub B6_ub B8_ub B9_ub B10_ub
    
   
    % 11th inequality constraint: 
    % biodiversity groups must experience 10% gain
%     if bio_constraint
%         num_groups = length(cnst_target);   
%         A_p = sparse(num_groups, num_options);
%         A_u = sparse(num_groups, num_farmers);
%         A_x = [];
%         for j = 1:num_farmers
%             A_x = [A_x, sparse(diag(cnst_data(j,:) - cnst_target'))];
%         end
%         A  = [A_p, A_u, A_x];
%         Bl = zeros(num_groups,1);
%         Bu = inf(num_groups,1);
%         if ~isempty(A),cplex.addRows(Bl, A, Bu); end
%         clear A A_p A_x Bl Bu    
%     end
    if any(cnst_target)
        num_groups = length(cnst_target); 
        cnst_data = reshape(cnst_data, length(cnst_target), []); % Reshape so have blocks of option to species group outcomes for each cell
        A_p = sparse(num_groups, num_options);
        A_u = sparse(num_groups, num_farmers);
        A_x = cnst_data;
        A  = [A_p, A_u, A_x];
        Bl = cnst_target;
        Bu = inf(num_groups,1);
        if ~isempty(A),cplex.addRows(Bl, A, Bu); end
        clear A A_p A_x Bl Bu    
    end
    
    
    % 2.1 Warm start
    % --------------
    if ~isempty(warm_start_prices)
        sln = [warm_start_prices warm_start_uptake];
        idx = [1:num_options num_options+num_farmers+1:num_options+num_farmers+num_options*num_farmers];
        idx = idx - 1;            
        filename = [cplex_options.logs 'warmstart.mst'];
        probname = 'elms_lp';
        fcn_write_warmstart(sln', idx', filename, probname);
        cplex.readMipStart(filename);
    end
    
    
    % 3. SOLVE
    % ========
    tic
    cplex.solve();
    toc
    
    
    % 4. OUTPUT
    % =========
    x        = cplex.Solution.x;
    prices   = x(1:num_options)';
    uptake   = round(reshape(x((num_options + num_farmers + 1):end), num_options, num_farmers)');
    fval     = cplex.Solution.objval;
    exitflag = cplex.Solution.status;
    exitmsg  = cplex.Solution.statusstring;
    
end



    