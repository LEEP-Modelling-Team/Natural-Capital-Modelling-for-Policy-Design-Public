% paymech_fr_activity_pctl_rnd.m
% ===============================
% Establishes prices at percentile of cost distribution and selects farmers
% at random from set who volunteer into an option at those prices. Repeats
% Niter times to assess probability of uptake across all farmers.

% 1. Initialise
% -------------
clear
rng(23112010)

% Model
% -----
payment_mechanism = 'fr_act_pctl_rnd';
unscaled_budget   = 1e9;
urban_pct_limit   = 0.5;
bio_constraint    = 0;
bio_as_prices     = false;  % only set to true if have a biodiversity const
byparcel          = true;
carbon_price_string = 'non_trade_central';
drop_vars = {'habitat_non_use', 'biodiversity'};
budget_str = [num2str(round(unscaled_budget/1e9)) 'bill'];
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];

pctl  = 50; % median prices
Niter = 1000;

% Markup
% ------
markup = 1.15;

% Paths to Data & Cplex Working Dir
% ---------------------------------
data_folder  = 'D:\myGitHub\defra-elms\Data\';
data_path = [data_folder, 'elm_data_', carbon_price_string, '.mat'];


% 2. Prepare data
% ---------------
data_year = 1;    
sample_size = 'no';  % all data
[b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);
num_prices  = length(price_vars);
num_options = size(b,2);
num_farmers = size(b,1);


% 3. Multiple Choice Knapsack Optimisation
% ----------------------------------------

% Calculate cost per hectare of each option
% -----------------------------------------
c(c==0) = inf;
c_perha = c./q;
c_perha(isinf(c_perha)) = NaN;
c_perha(isnan(c_perha)) = NaN;

% Prices at Percentile
% --------------------
prices = prctile(c_perha, pctl, 1);

% Uptake and Benefits at those prices
% -----------------------------------
uptake   = myfun_uptake(prices, q, c, elm_options);    
benefits = sum(b.*uptake, 2);
spend    = sum(prices.*q.*uptake, 2);

% Number of Iterations of Random Uptake Selection
% -----------------------------------------------
mc_spend  = zeros(Niter, 1);
mc_fval   = zeros(Niter, 1);
mc_uptake = zeros(num_farmers, 1);
mc_cnst   = zeros(Niter, 1);

for i = 1:Niter
    
    % Random first come, first served
    % -------------------------------
	sortidx = randperm(num_farmers);
	inbudget_ind = 1 - (cumsum(spend(sortidx)) >= budget);
	inbudget_ind(sortidx) = inbudget_ind;  % undo sort to find best in budget
 
    % Select most benefits in budget
    % ------------------------------
    uptake_i    = uptake .* inbudget_ind;
    mc_uptake   = mc_uptake + sum(uptake_i,2);        
    mc_spend(i) = sum(prices.*q.*uptake_i, 'all');
    mc_fval(i)  = sum(b.*uptake_i, 'all'); 
        
    % Check if meets constraint
    % -------------------------
    num_spgrp = length(cnst_target);
    spgrp_chg = zeros(num_spgrp,1);
    for k = 1:num_spgrp
        spgrp_chg(k) = sum(uptake_i.*squeeze(cnst_data(k,:,:))', 'all');        
    end
    if ~any(spgrp_chg < cnst_target)
        mc_cnst(i) = 1;
    end
        
end

% Uptake as a probability
% -----------------------
rnd_uptake_ind = mc_uptake/Niter;

% Constraint satisfaction as a probability
% ----------------------------------------
cnst_prob = sum(mc_cnst)/Niter;

% Process result
% --------------
uptake        = uptake .* rnd_uptake_ind;
option_nums   = (1:8)';
option_choice = (uptake>0) * option_nums;
benefits      = sum(b.*uptake, 2);
costs         = nansum(c.*uptake, 2);
farm_payment  = sum(prices.*q.*uptake, 2);
hectares      = sum(table2array(hectares).*uptake, 2);

% 4. Save Solution
% ----------------
solution.prices        = prices;
solution.fval          = sum(benefits);
solution.spend         = sum(farm_payment);
solution.cnst_prob     = cnst_prob;
solution.uptake        = uptake;
solution.uptake_ind    = rnd_uptake_ind;
solution.option_choice = option_choice;
solution.new2kid       = new2kid(rnd_uptake_ind>0);
solution.hectares      = hectares;
solution.farm_costs    = costs;
solution.farm_benefits = benefits;
solution.farm_payment  = farm_payment;


if bio_constraint > 0    
    save(['solution_' payment_mechanism  '_' budget_str '_' biocnst_str '.mat'], 'solution'); 
else
    save(['solution_' payment_mechanism  '_' budget_str '.mat'], 'solution');     
end




