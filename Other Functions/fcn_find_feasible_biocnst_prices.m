% fcn_find_feasible_biocnst_prices
% ================================

% Purpose
% -------
%  Search price space to find a set of feasible prices that do not break
%  the constraints (budget & biodiverisity targets)

function [prices_feas, spend_feas] = fcn_find_feasible_biocnst_prices(budget, c, q, elm_options, prices_uncnst, prices_lb, prices_ub, cnst_data, cnst_target)
        
    % Constants
    % ---------
    num_spgrp = length(cnst_target);
        
    % Check for constraint violation
    % ------------------------------
    uptake = myfun_uptake(prices_uncnst, q, c, elm_options);  
    spgrp_chg = zeros(num_spgrp,1);
    for k = 1:num_spgrp
        spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
    end
    % [spgrp_chg cnst_target (spgrp_chg-cnst_target)>0]    
    
    % Search for Feasible Prices
    % --------------------------
    if all((spgrp_chg-cnst_target) > 0) 
        prices_feas = prices_uncnst;
    else        
        % Scale so prices same order of magnitude for maximisation
        % --------------------------------------------------------
        prices_uncnst(prices_uncnst==0) = 1;
        ord_mag = 10.^floor(log(abs(prices_uncnst))./log(10));
        ord_mag(ord_mag==0) = 1;
        prices_uncnst = prices_uncnst./ord_mag;

        q_scl = q.*ord_mag;
        prices_lb = prices_lb ./ ord_mag';
        prices_lb(isnan(prices_lb)) = 0;
        prices_ub = prices_ub ./ ord_mag'; 

        % GA search for biodiversity constrained feasible prices
        % ------------------------------------------------------
        % Objective function is to minimise spend subject to not violating
        % the biodiversity constraints with a stopping function that at the
        % end of each generation checks to see if the best objective value
        % has dropped below the budget.
        benefitfunc    = @(p) myfun_spend(p, q_scl, c, elm_options);
        constraintfunc = @(p) mycon_Biod(p, q_scl, c, elm_options, cnst_data, cnst_target);   
        stopfunc       = @(options, state, flag) mystop_budget(options, state, flag, q_scl, c, elm_options, budget, cnst_data, cnst_target);   

        tic
        options = optimoptions('ga', ...
                               'OutputFcn', stopfunc, ...
                               'Display', 'iter');
        [prices_feas, spend_feas] = ga(benefitfunc, length(prices_uncnst),[],[],[],[],prices_lb,prices_ub,constraintfunc,options);
        toc  
        
        % Does the solution break the budget contraint?
        % ---------------------------------------------
        if spend_feas > budget
            warning('Unable to find prices that can deliver biodiversity constraint in budget');
        else
            prices_feas = prices_feas.*ord_mag;
        end
    end
    
 end