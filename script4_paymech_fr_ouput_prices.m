% paymech_opt_output_prices.m
% ===============================
%  Search for optimum flat rate prices for outcomes to deliver optimal 
%  benefits for budget.

% 1. Initialise
% -------------
clear
rng(23112010)

% Reuse Previous Prices and/or Solutions
% --------------------------------------
use_old_soln       = false;
use_old_prices_sml = false;

% Model
% -----
payment_mechanism = 'fr_es';
unscaled_budget   = 1e9;
urban_pct_limit   = 0.5;
byparcel          = true;
sample_size       = 5000;
carbon_price_string = 'non_trade_central';
drop_vars   = {'habitat_non_use', 'biodiversity'};
budget_str  = [num2str(round(unscaled_budget/1e9)) 'bill' ];
if sample_size > 1000            
	sample_str = [num2str(round(sample_size/1000)) 'k_sample'];
else
    sample_str = [num2str(round(sample_size)) '_sample'];
end

% Markup
% ------
markup = 1.15;

% Paths to Data, Solutions & Cplex Working Dir
% --------------------------------------------
cplex_folder  = 'D:\myGitHub\defra-elms\Cplex\'; % Path to use for cplex working files
data_folder   = 'D:\myGitHub\defra-elms\Data\'; % Path to data
soln_folder   = 'D:\myGitHub\Natural-Capital-Modelling-for-Policy-Design\'; % Path to stored old solution
prices_folder = 'D:\myGitHub\Natural-Capital-Modelling-for-Policy-Design\'; % Path to stored price matrix
data_path = [data_folder, 'elm_data_', carbon_price_string, '.mat'];



% 2. Unconstrained Problem
% ------------------------
%   The full optimisation problem is too large to be solved directly as a
%   MIP. We employ a number of strategies to identify a feasible price
%   vector to act as a warm start for the MIP and to find upper and lower
%   bounds for prices to limit search.
bio_constraint     = 0;      % Just budget constrained
bio_as_prices      = false;  % No biodiversity prices
use_old_soln       = false;
use_old_prices_sml = false;
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];


% if bio_constraint == 0
% 
%     % Prices matrix file name 
%     % -----------------------
%     pricesfile_name = [prices_folder 'prices_' payment_mechanism '_' budget_str '_' biocnst_str '_' sample_str '.mat'];
% 
%     if use_old_soln
% 
%         % 2.1 Use Previously-Estimated Solution
%         % -------------------------------------
%         load([soln_folder 'solution_' payment_mechanism '_' budget_str '.mat'], 'solution'); 
%         prices_ga = solution.prices;
%         uptake_ga = solution.uptake;
% 
%         load(pricesfile_name, "prices_sml")
% 
%     else    
% 
%         % 2.2 Prices for Small Samples
%         % ----------------------------
%         if use_old_prices_sml
% 
%             % Load previously-estimated matrix of prices from small samples
%             % -------------------------------------------------------------
%             load(pricesfile_name, "prices_sml")    
%         else
%             % Estimate prices from small samples
%             % ----------------------------------
%             Niter = 10;
%             [prices_sml, benefits_sml] = small_sample_prices(Niter, sample_size, unscaled_budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);       
%             save(pricesfile_name, "prices_sml", "benefits_sml")
%         end
% 
%         % 2.3 Full Sample: GA Solve
%         % -------------------------
%         prices = [];
%         [prices_ga, benefit_ga] = full_sample_ga(prices, prices_sml, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
% 
%     end
% 
%     % 2.4 Full Sample: MIP Solve
%     % --------------------------
%     solution = full_sample_mip(prices_ga, prices_sml, unscaled_budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
% 
%     save([soln_folder 'solution_' payment_mechanism '_' budget_str '.mat'], 'solution');  
% 
% end



% 3. Biodiversity Constrained Problem
% -----------------------------------
bio_constraint = 0.15;    % 15% uplift in biodiversity
bio_as_prices  = false;   % No biodiversity prices
use_old_soln   = false;  
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];

% if bio_constraint > 0 && ~bio_as_prices
% 
%     % Prices matrix file name 
%     % -----------------------
%     pricesfile_name = [prices_folder 'prices_' payment_mechanism '_' budget_str '_0pct_' sample_str '.mat'];
% 
%     if use_old_soln
% 
%         % 3.1 Use Previously-Estimated Solution
%         % -------------------------------------        
%         load([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
%         prices_ga = solution.prices;
%         uptake_ga = solution.uptake;
% 
%         load(pricesfile_name, "prices_sml");
% 
%         % 3.2 Full Sample: MIP Solve
%         % --------------------------
%         solution = full_sample_mip(prices_ga, prices_sml, budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
%         save([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution');         
% 
%     else    
% 
%         % Load data
%         % ---------
%         data_year = 1;    
%         sample_size = 'no';  % all data
%         [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
% 
%         % 3.3 Unconstrained Solution
%         % --------------------------
% 
%         % Load solution
%         % -------------
%         load([soln_folder 'solution_' payment_mechanism '_' budget_str '.mat'], 'solution');        
%         prices_uncnst = solution.prices;
% 
%         % Load prices used to find solution
%         % ---------------------------------
%         load([prices_folder 'prices_' payment_mechanism '_' budget_str '_0pct_' sample_str '.mat'], 'prices_sml')        
% 
% 
%         % Check for constraint violation
%         % ------------------------------
%         uptake = myfun_uptake(prices_uncnst, q, c, elm_options);
%         num_spgrp = length(cnst_target);
%         spgrp_chg = zeros(num_spgrp,1);
%         for k = 1:num_spgrp
%             spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
%         end
% 
%         if any((spgrp_chg-cnst_target) < 0) % Unconstrained Solution fails Biodiversity Constrained Problem
% 
%             % 3.4 Check Problem can be Solved
%             % -------------------------------
%             %  Can any prices meet the biodiversity constraints in the
%             %  budget?
% 
%             % Price limits for search
%             % -----------------------
%             prices_max = fcn_find_max_prices(q, c, budget, lu_data, elm_options);
%             prices_lb = zeros(size(prices_max)); 
%             prices_ub = max(prices_sml);
% 
%             % Reduce to relevant cells
%             % ------------------------
%             %  Remove cells where no possible price would enduce participation
%             excluded_cells = fcn_find_unusable_cells(q, c, budget, elm_options, prices_max, 1.5*[prices_uncnst; prices_ub]');
%             b(excluded_cells, :)          = [];
%             c(excluded_cells, :)          = [];
%             q(excluded_cells, :, :)       = [];
%             lu_data(excluded_cells, :)    = []; 
%             cnst_data(:,:,excluded_cells) = [];        
% 
%             % Search for Feasible Prices
%             % --------------------------
%             [prices_feas, spend_score] = fcn_find_feasible_biocnst_prices(budget, c, q, elm_options, prices_uncnst, prices_lb, prices_ub', cnst_data, cnst_target);
%             if spend_score > budget
%                 error('No feasible prices to deliver biodiversity constraint in budget');
%             end
% 
%             % 3.5 Full Sample: GA Solve
%             % -------------------------
%             [prices_ga, benefit_ga] = full_sample_ga([prices_feas; prices_uncnst], prices_sml, budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
% 
%             % 3.6 Full Sample: MIP Solve
%             % --------------------------
%             solution = full_sample_mip([prices_ga; prices_feas; prices_uncnst], prices_sml, budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
%             save([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
% 
%         else
% 
%             % 3.7 Unconstrained Solution solves Biodiversity Constrained Problem
%             % ------------------------------------------------------------------    
%             save([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
%             save(['solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution');                    
%             copyfile([soln_folder  'prices_' payment_mechanism '_' budget_str '_0pct_' sample_str '.mat'], ...
%                      [soln_folder  'prices_' payment_mechanism '_' budget_str '_' biocnst_str '_' sample_str '.mat']);   
% 
%         end
% 
%     end
% 
% end



% 4. Biodiversity Constrained Problem with Biodiversity Prices
% ------------------------------------------------------------
bio_constraint = 0.15;    % 15% uplift in biodiversity
bio_as_prices  = true;    % Use biodiversity prices
use_old_soln   = false;  
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];

if bio_constraint > 0 && bio_as_prices

    % Load data
    % ---------
    data_year = 1;    
    sample_size = 'no';  % all data
    [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
    num_spgrp = length(cnst_target);

    if use_old_soln

        % 4.1 Use Previously-Estimated Solution
        % -------------------------------------        
        load([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '_pbio.mat'], 'solution'); 
        prices_mip = solution.prices;     
        
        prices_outcomes = solution.prices(1:end-num_spgrp);
        prices_biod     = solution.prices(end-num_spgrp+1:end);

        prices_ub  = [prices_outcomes, ones(1, num_spgrp)*100;
                      prices_outcomes, prices_biod*3];

        % 4.2 Full Sample: MIP Solve
        % --------------------------
        solution = full_sample_mip(prices_mip, prices_ub, budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
        save([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '_pbio.mat'], 'solution');

    else    

        % 4.3 Solution without Biodiversity P+rices
        % ----------------------------------------
        load([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution');        
        prices_feas = [solution.prices, zeros(1, num_spgrp)];                   
        prices_ub   = [solution.prices, ones(1, num_spgrp)*100;
                       solution.prices, ones(1, num_spgrp)*100];

        % 4.4 Full Sample: GA Solve
        % -------------------------
        [prices_ga, benefit_ga] = full_sample_ga(prices_feas, prices_ub, budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);

        % 4.5 Full Sample: MIP Solve
        % --------------------------
        solution = full_sample_mip([prices_ga; prices_feas], prices_ub, budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel);
        save([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '_pbio.mat'], 'solution');

    end
    
end



% Function: small_sample_prices
% -----------------------------
%  Takes repeated small samples from data and solves for prices using a ga
%  algorithm and then a MIP.

function [prices, benefits] = small_sample_prices(Niter, sample_size, unscaled_budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year)
    
    % 1. Allocate Small Sample Solution Matrices
    % ------------------------------------------
    prices   = [];
    benefits = []; 

    for iter = 1:Niter

        fprintf('Iteration: %d of %d\n', iter, Niter);
        fprintf('------------------\n');

        % 2. Load new sample of data
        % --------------------------
        data_year = 1;    % year in which scheme run 
        [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
        num_farmers = size(q, 1);
        num_prices  = size(q, 2);
        num_options = size(q, 3);

        % 3. Maximum possible prices
        % --------------------------
        prices_max = fcn_find_max_prices(q, c, budget, lu_data, elm_options);
        prices_min = zeros(size(prices_max));   

        % 4. Scale quantities
        % -------------------
        % Scale so prices same order of magnitude for maximisation
        ord_mag = 10.^floor(log(abs(prices_max'))./log(10));
        q_scl = zeros(size(q));
        for i = 1:length(elm_options)
            q_scl(:, :, i) = q(:, :, i) .* ord_mag;
        end
        prices_min = prices_min ./ ord_mag';
        prices_max = prices_max ./ ord_mag'; 

        % 5. Reduce to relevant cells
        % ---------------------------
        %  Remove cells where no possible price would enduce participation
        excluded_cells = fcn_find_unusable_cells(q_scl, c, budget, elm_options, prices_max);
        b(excluded_cells, :)          = [];
        c(excluded_cells, :)          = [];
        q_scl(excluded_cells, :, :)   = [];
        lu_data(excluded_cells, :)    = []; 
        cnst_data(:,:,excluded_cells) = [];

        % 6. Solution: ga algorithm
        % -------------------------
        tic
        benefitfunc    = @(p) myfun_ES(p, q_scl, c, b, elm_options);
        if ~bio_constraint
            constraintfunc = @(p) mycon_budget(p, q_scl, c, budget, elm_options);
        else
            constraintfunc = @(p) mycon_budget_bio(p, q_scl, c, budget, elm_options, cnst_data, cnst_target);
        end
        options = optimoptions('ga', ...
                               'Display', 'iter');
        [prices_ga, benefit_ga] = ga(benefitfunc, num_prices,[],[],[],[],prices_min,prices_max,constraintfunc,options);
        toc
        uptake_ga = myfun_uptake(prices_ga, q_scl, c, elm_options)';
        uptake_ga = uptake_ga(:);

        % 7. Scale prices to ensure do not violate budget 
        % -----------------------------------------------
        %  Want to ensure that this is a feasible solution to act as a
        %  warmstart in the MIP
        spend = myfun_spend(prices_ga, q_scl, c, elm_options);
        if spend > budget
            budgetgapfun = @(scale)((budget-2e5-myfun_spend(prices_ga*scale, q_scl, c, elm_options))^2);
            budgetcon    = @(scale)(myfun_spend(prices_ga*scale, q_scl, c, elm_options)-budget);
            budgetconfun = @(scale)deal(budgetcon(scale), []);
            options = optimoptions('fmincon', ...
                                   'ConstraintTolerance', 0, ...
                                   'Display', 'none');
            scale = fmincon(budgetgapfun,budget/spend, [], [], [], [], 0, 10*budget/spend, budgetconfun, options);
            % Rescale prices and uptake
            prices_ga = prices_ga * scale;
            uptake_ga = myfun_uptake(prices_ga, q_scl, c, elm_options)';
            uptake_ga = uptake_ga(:);
        end

        % 8. Solution: MIP
        % ----------------
        cplex_options.time = 800;
        cplex_options.logs = cplex_folder;    
        % Restrict price search space
        prices_min = prices_ga' * 0.5;
        prices_max = prices_ga' * 1.5;
        [prices_sml, uptake_sml, fval, exitflag, exitmsg] = MIP_fr_out(b, c, q_scl, budget, lu_data, prices_ga', uptake_ga, prices_min, prices_max, cnst_data, cnst_target, byparcel, cplex_options);

        % 9. Solution: Store 
        % ------------------
        if exitflag == 101 || exitflag == 102 || exitflag == 107
            prices   = [prices; prices_sml .* ord_mag];
            benefits = [benefits; fval];
        end

    end   

    % 10. Scale prices_sml to clear full budget
    % -----------------------------------------
    sample_size = 'no';  % all data
    data_year = 1;    
    [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
    
    for i = 1:height(prices)

        % Find scale of prices to fit full problem within budget
        % ------------------------------------------------------
        spend = myfun_spend(prices(i,:), q, c, elm_options);
        budgetgapfun = @(scale)((budget-2e6-myfun_spend(prices(i,:)*scale, q, c, elm_options))^2);
        budgetcon    = @(scale)(myfun_spend(prices(i,:)*scale, q, c, elm_options)-budget);
        budgetconfun = @(scale)deal(budgetcon(scale), []);
        options = optimoptions('fmincon', ...
                               'ConstraintTolerance', 0, ...
                               'Display', 'none');
        scale = fmincon(budgetgapfun,budget/spend, [], [], [], [], 0, 10*budget/spend, budgetconfun, options);
        % myfun_spend(prices(i,:)*scale, q, c, elm_options)

        % Update prices and benefits to full problem solution
        % ---------------------------------------------------
        prices(i,:) = prices(i,:) * scale;
        benefits(i) = -myfun_ES(prices(i,:), q, c, b, elm_options);

    end

end


% Function: full_sample_ga
% ------------------------
%  Solves for the full sample using a ga algorithm initiated from a matrix
%  of prices estimated from a set of small samples.
function [prices_ga, benefit_ga] = full_sample_ga(prices, prices_sml, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year)

    % 1. Load data
    % ------------
    sample_size = 'no';  % all data
    data_year = 1;    
    [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
    num_farmers = size(q, 1);
    num_prices  = size(q, 2);
    num_options = size(q, 3);
    
    % 2. Price bounds from sample searches
    % ------------------------------------    
    prices_max = fcn_find_max_prices(q, c, budget, lu_data, elm_options);
    prices_lb  = zeros(size(prices_max)); 
    prices_ub  = max(prices_sml)';

    % 3. Scale quantities
    % -------------------
    % Scale so prices same order of magnitude for maximisation
    ord_mag = 10.^floor(log(abs(prices_ub'))./log(10));
    ord_mag(ord_mag==0) = 1;
    q_scl = zeros(size(q));
    for i = 1:length(elm_options)
        q_scl(:, :, i) = q(:, :, i).*ord_mag;
    end
    prices_lb  = prices_lb  ./ ord_mag';
    prices_ub  = prices_ub  ./ ord_mag'; 
    prices_max = prices_max ./ ord_mag'; 
    prices     = prices     ./ ord_mag;
    prices_sml = prices_sml ./ ord_mag;
    
    prices_lb(isinf(prices_lb)) = 0;
    
    % 5. Reduce to relevant cells
    % ---------------------------
    %  Remove cells where no possible price would enduce participation
    excluded_cells = fcn_find_unusable_cells(q_scl, c, budget, elm_options, prices_max, 1.5*[prices', prices_ub]);
    b(excluded_cells, :)          = [];
    c(excluded_cells, :)          = [];
    q_scl(excluded_cells, :, :)   = [];
    lu_data(excluded_cells, :)    = []; 
    cnst_data(:,:,excluded_cells) = [];
    
    % 6. Solve with ga algorithm
    % --------------------------
    % Benefit and Constraint Functions for ga
    benefitfunc    = @(p) myfun_ES(p, q_scl, c, b, elm_options);
    if ~bio_constraint
        constraintfunc = @(p) mycon_budget(p, q_scl, c, budget, elm_options);
    else
        constraintfunc = @(p) mycon_budget_bio(p, q_scl, c, budget, elm_options, cnst_data, cnst_target);
    end
    % sml_prices as starting population
    options = optimoptions('ga', ...
                           'InitialPopulationMatrix', [prices; prices_sml], ... 
                           'PopulationSize', height([prices; prices_sml]), ...
                           'Display', 'iter');
    [prices_ga, benefit_ga] = ga(benefitfunc, num_prices,[],[],[],[], prices_lb, prices_ub*1.5, constraintfunc,options);
        
    % 7. Solution
    % -----------
    prices_ga  = prices_ga  .* ord_mag;

end


% Function: full_sample_mip
% ------------------------
%  Solves for the full sample using a mip algorithm initiated from a matrix
%  of prices estimated from a set of small samples and a warm start 
%  estimated using a ga algorithm.
function solution = full_sample_mip(prices, prices_sml, unscaled_budget, data_path, cplex_folder, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year)

    % 1. Load data
    % ------------
    sample_size = 'no';  % all data
    data_year = 1;    
    [b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
    num_farmers = size(q, 1);
    num_prices  = size(q, 2);
    num_options = size(q, 3);
    
    % 2. Price bounds from sample searches
    % ------------------------------------  
    prices_max = fcn_find_max_prices(q, c, budget, lu_data, elm_options);
    prices_lb  = zeros(size(prices_max));
    prices_ub  = max(prices_sml)';
        
    % 3. Scale quantities
    % -------------------
    % Scale so prices same order of magnitude for maximisation
    ord_mag = 10.^floor(log(abs(prices_ub'))./log(10));
    ord_mag(ord_mag==0) = 1;
    q_scl = zeros(size(q));
    for i = 1:length(elm_options)
        q_scl(:, :, i) = q(:, :, i).*ord_mag;
    end
    prices_lb  = prices_lb  ./ ord_mag';
    prices_ub  = prices_ub  ./ ord_mag'; 
    prices_max = prices_max ./ ord_mag'; 
    prices     = prices     ./ ord_mag;
    prices_sml = prices_sml ./ ord_mag;
   
    prices_lb(isinf(prices_lb)) = 0;
    
    % 4. Reduce to relevant cells
    % ---------------------------
    %  Remove cells where range of plausible prices would enduce participation
    excluded_cells = fcn_find_unusable_cells(q_scl, c, budget, elm_options, prices_max, 1.5*[prices', prices_ub]);
    b(excluded_cells, :)          = [];
    c(excluded_cells, :)          = [];
    q_scl(excluded_cells, :, :)   = [];
    lu_data(excluded_cells, :)    = []; 
    cnst_data(:,:,excluded_cells) = [];
    
    % 5. Uptake for Warmstart
    % -----------------------
    prices_start = [];
    uptake_start = [];
    for i = 1:height(prices)
        prices_start = [prices_start  prices(i,:)'];
        uptake_start = [uptake_start reshape(myfun_uptake(prices(i,:), q_scl, c, elm_options)', [], 1)]; 
    end
    for i = 1:height(prices_sml)
        prices_start = [prices_start prices_sml(i,:)'];
        uptake_start = [uptake_start reshape(myfun_uptake(prices_sml(i,:), q_scl, c, elm_options)', [], 1)]; 
    end
    
    % 6. Solve MIP
    % ------------
    cplex_options.time = 15000;
    cplex_options.logs = cplex_folder;    
    [prices, uptake_sml, fval, exitflag, exitmsg] = MIP_fr_out(b, c, q_scl, budget, lu_data, prices_start, uptake_start, prices_lb*0.5, prices_ub*1.5, cnst_data, cnst_target, byparcel, cplex_options);
        
    % 7. Process result
    % -----------------
    uptake_ind_sml    = (sum(uptake_sml,2) > 0);
    option_nums       = (1:8)';
    option_choice_sml = (uptake_sml * option_nums);
    benefits_sml      = sum(b.*uptake_sml, 2);
    costs_sml         = sum(c.*uptake_sml, 2);
    payments          = zeros(size(uptake_sml));
    for i = 1:num_options
        payments(:,i)  = prices * q_scl(:, :, i)';
    end
    farm_payment_sml = sum(payments.*uptake_sml, 2);
    
    % 8. Re-expand to full cell list
    % ------------------------------
    uptake        = zeros(num_farmers, num_options);
    uptake_ind    = zeros(num_farmers, 1);
    option_choice = zeros(num_farmers, 1);
    benefits      = zeros(num_farmers, 1);
    costs         = zeros(num_farmers, 1);
    farm_payment  = zeros(num_farmers, 1);
    
    sample_idx                = find(1-excluded_cells);
    uptake(sample_idx,:)      = uptake_sml;
    uptake_ind(sample_idx)    = uptake_ind_sml;
    uptake_ind                = logical(uptake_ind);
    option_choice(sample_idx) = option_choice_sml;
    benefits(sample_idx)      = benefits_sml;
    costs(sample_idx)         = costs_sml;
    farm_payment(sample_idx)  = farm_payment_sml;
    
    % 9. Solution
    % -----------
    solution.prices        = prices .* ord_mag;
    solution.fval          = sum(benefits);
    solution.spend         = sum(farm_payment);
    solution.uptake        = uptake;
    solution.uptake_ind    = uptake_ind;
    solution.option_choice = option_choice;
    solution.new2kid       = new2kid(uptake_ind);
    solution.hectares      = full(sum(table2array(hectares).*uptake,2));
    solution.farm_costs    = costs;
    solution.farm_benefits = benefits;
    solution.farm_payment  = farm_payment;

end