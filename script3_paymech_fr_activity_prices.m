% paymech_fr_activity_prices.m
% =============================
%  Search for optimum flat rate prices for activities to deliver optimal 
%  benefits for budget.

% 1. Initialise
% -------------
clear
rng(23112010)

% Reuse Previous Prices and/or Solutions
% --------------------------------------
use_old_soln = true;
use_old_sml_prices = true;

% Model
% -----
payment_mechanism = 'fr_act';
unscaled_budget   = 1e9;
urban_pct_limit   = 0.5;
bio_constraint    = 0.15;
bio_as_prices     = true;
byparcel          = true;
carbon_price_string = 'non_trade_central';
drop_vars = {'habitat_non_use', 'biodiversity'};
budget_str  = [num2str(round(unscaled_budget/1e9)) 'bill'];
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];

% Markup
% ------
markup = 1.15;

% Paths to Data, Solutions & Cplex Working Dir
% --------------------------------------------
cplex_folder = 'D:\myGitHub\defra-elms\Cplex\'; % Path to use for cplex working files
data_folder  = 'D:\myGitHub\defra-elms\Data\'; % Path to data
soln_folder  = 'D:\myGitHub\Natural-Capital-Modelling-for-Policy-Design\Data\'; % Path to stored old solution
data_path = [data_folder, 'elm_data_', carbon_price_string, '.mat'];


% 2. Prepare data
% ---------------

% 2.1 Load data
% -------------
data_year = 1;    
sample_size = 'no';  % all data
[b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
num_prices  = length(price_vars);
num_options = size(b,2);
num_farmers = size(b,1);

% 2.2 Reduce problem size
% -----------------------
%  Reduce size of the problem: search for prices that would exhaust the
%  entire budget for each option, calculate uptake and remove subset of
%  cells that are never selected under any option when the budget is
%  exhausted
if ~bio_as_prices
    prices_lb = zeros(length(elm_options), 1);
    prices_ub = zeros(length(elm_options), 1);
    excluded_cells = ones(num_farmers,1);    
    for i = 1:length(elm_options)
        [price_pts, sortidx] = sort(c(:,i)./q(:,i), 'ascend', 'MissingPlacement', 'last');
        price_pts(isnan(price_pts)) = inf;
        ind_over_budget = (cumsum(q(sortidx,i)) .* price_pts) >= budget;
        cell_over_budget(sortidx) = ind_over_budget;
        excluded_cells = excluded_cells .* cell_over_budget';
        prices_ub(i)   = price_pts(find(ind_over_budget, 1)-1);
    end
else           
    prices_ub = fcn_find_max_prices(q, c, budget, lu_data, elm_options);
    prices_lb = zeros(size(prices_ub));
    excluded_cells = fcn_find_unusable_cells(q, c, budget, elm_options, prices_ub);
end

% 2.3 Reduce to relevant cells
% ----------------------------
%  Remove cells where no possible price would enduce participation
excluded_cells = logical(excluded_cells);
b(excluded_cells, :) = [];
c(excluded_cells, :) = [];
if ~bio_as_prices
    q(excluded_cells, :) = [];
else
    q(excluded_cells, :, :) = [];
end
lu_data(excluded_cells, :) = []; 
cnst_data(:,:,excluded_cells) = [];

% 2.4 Search for feasible prices with biodiversity constraint
% -----------------------------------------------------------
warm_start_prices = [];
warm_start_uptake = [];

if use_old_soln 

    % Load Solution for this Payment Mechanism
    % ----------------------------------------
    if bio_constraint > 0
        if ~bio_as_prices
            load([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
        else
            load([soln_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '_pbio.mat'], 'solution'); 
        end
    else
        load([soln_folder 'solution_' payment_mechanism '_' budget_str '.mat'], 'solution');     
    end    

    % Remove excluded farms
    uptake = solution.uptake;
    uptake(excluded_cells, :) = []; 
            
    warm_start_prices = solution.prices;
    warm_start_uptake = uptake';
    warm_start_uptake = warm_start_uptake(:)';    

else
    if bio_constraint > 0
            
        if bio_as_prices 
            % Solution prices from biodiversity constrained model (without biodiversity prices) 
            % (Assumes already run the constrained model if trying to add biodiversity prices)
            load([data_folder 'solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution');        
            prices_feas = solution.prices;          
            prices_feas = [prices_feas, zeros(1, length(cnst_target))];        
                    
            % Check for constraint violation
            % ------------------------------
            % uptake = myfun_uptake(prices_feas, q, c, elm_options);
            uptake = solution.uptake(~excluded_cells, :);
            num_spgrp = length(cnst_target);
            spgrp_chg = zeros(num_spgrp,1);
            for k = 1:num_spgrp
                spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
            end
            
            % Search for Feasible Prices 
            % --------------------------        
            if any((spgrp_chg-cnst_target) < 0) 
                % Constrained prices should be feasible (they meet biodiversity
                % and budget constraints), but they may fail as Cplex breaks 
                % ties in the profitably of options in a cell in favour of 
                % constraint satisfaction, but that logic is difficult to 
                % implement in simple matlab code for constraint checks.            
                benefitfunc    = @(p) myfun_ES(p, q, c, b, elm_options);
                constraintfunc = @(p) mycon_budget_bio(p, q, c, budget, elm_options, cnst_data, cnst_target);
                % Base initial population on constrained prices. 
                popSize = 50; % population size
                initialpop = lhsdesign(popSize, num_prices); % generate LHS sample
                initialpop = prices_feas*0.9 + initialpop .* (prices_feas*1.1 - prices_feas*0.9);
                initialpop = [prices_feas; initialpop];                      
                tic
                options = optimoptions('ga', ...
                                       'InitialPopulation', initialpop, ...
                                       'Display', 'iter');
                [prices_feas, benefit_feas] = ga(benefitfunc, length(prices_feas),[],[],[],[],prices_lb,prices_ub,constraintfunc,options);
                toc  
                if isempty(prices_feas)
                    error('No feasible prices to deliver biodiversity constraint in budget');
                end
            end
            
            warm_start_prices = prices_feas;
            warm_start_uptake = uptake';
            warm_start_uptake = warm_start_uptake(:)';
            
        else
            % Solution prices from unconstrained model
            % ----------------------------------------
            load([data_folder 'solution_' payment_mechanism '_' budget_str '.mat'], 'solution');        
            prices_uncnst = solution.prices;
    
            % Check for constraint violation
            % ------------------------------
            uptake = myfun_uptake(prices_uncnst, q, c, elm_options);
            num_spgrp = length(cnst_target);
            spgrp_chg = zeros(num_spgrp,1);
            for k = 1:num_spgrp
                spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
            end
    
            % Search for Feasible Prices
            % --------------------------
            if any((spgrp_chg-cnst_target) < 0) 
    
                prices_feas = fcn_find_feasible_biocnst_prices(budget, c, q, elm_options, prices_uncnst, prices_lb, prices_ub, cnst_data, cnst_target);
    
                if isempty(prices_feas)
                    error('No feasible prices to deliver biodiversity constraint in budget');
                else
                    
                    % Check for constraint violation
                    % ------------------------------
                    uptake = myfun_uptake(prices_feas, q, c, elm_options);
                    num_spgrp = length(cnst_target);
                    spgrp_chg = zeros(num_spgrp,1);
                    for k = 1:num_spgrp
                        spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
                    end
                    [cnst_target spgrp_chg ((spgrp_chg-cnst_target) > 0)]                
                    
                    warm_start_prices = prices_feas;
                    warm_start_uptake = myfun_uptake(prices_feas, q, c, elm_options)';
                    warm_start_uptake = warm_start_uptake(:)';
                end
            else
                warm_start_prices = prices_feas;
                warm_start_uptake = uptake';
                warm_start_uptake = warm_start_uptake(:)';
            end
    
        end
    end
end


% 3. MIP for Global Optimal Prices
% --------------------------------
cplex_options.time = 29000;
cplex_options.logs = cplex_folder;

if ~bio_as_prices 
    [prices, uptake_sml, fval, exitflag, exitmsg] = MIP_fr_act(b, c, q, budget, lu_data, warm_start_prices, warm_start_uptake, prices_lb, prices_ub, cnst_data, cnst_target, cplex_options);
else
    [prices, uptake_sml, fval, exitflag, exitmsg] = MIP_fr_out(b, c, q, budget, lu_data, warm_start_prices, warm_start_uptake, prices_lb, prices_ub, cnst_data, cnst_target, byparcel, cplex_options);
end
        
% 3.1 Process result
% ------------------
uptake_ind_sml    = (sum(uptake_sml,2) > 0);
option_nums       = (1:8)';
option_choice_sml = (uptake_sml * option_nums);
benefits_sml      = sum(b.*uptake_sml, 2);
costs_sml         = sum(c.*uptake_sml, 2);
if ~bio_as_prices 
    farm_payment_sml  = sum(prices.*q.*uptake_sml, 2);
else
    payments = zeros(size(uptake_sml));
    for i = 1:num_options
        payments(:,i)  = prices * q(:, :, i)';
    end
    farm_payment_sml = sum(payments.*uptake_sml, 2);
end

% 3.2 Re-expand to full cell list
% -------------------------------
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


% 4. Save Solution
% ----------------
solution.prices        = prices;
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

if bio_constraint > 0
    if ~bio_as_prices
        save(['solution_' payment_mechanism '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
    else
        save(['solution_' payment_mechanism '_' budget_str '_' biocnst_str '_pbio.mat'], 'solution'); 
    end
else
    save(['solution_' payment_mechanism '_' budget_str '.mat'], 'solution');     
end
