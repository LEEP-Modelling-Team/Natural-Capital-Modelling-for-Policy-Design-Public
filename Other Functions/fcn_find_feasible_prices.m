% fcn_find_feasible_prices
% ========================

% Purpose
% -------
%  Search price space to find a set of feasible prices that do not break
%  the constraints (budget & biodiverisity targets)

function [prices_feasible, benefits_feasible] = fcn_find_feasible_prices(budget, benefits, costs, q, elm_options, prices_min, prices_max, cnst_data, cnst_target, Nfeas)

    % fprintf('\n  Rough Search of Parameter Space:');
    % fprintf('\n  ------------------------------- \n');
    N  = 15000;
    Nq = size(q,2);
    
    % Choose some test prices
    % -----------------------
    LatinHC     = lhsdesign(N, Nq);
    prices_test = prices_min + LatinHC .* (prices_max - prices_min);
    
    % Evaluate benefits and spend at test prices
    % ------------------------------------------
    benefits_test = nan(N,1);
    parfor i = 1:N
        benefits_test(i) = -myfun_ES(prices_test(i, :), q, costs, benefits, elm_options);
        if any(mycon_budget_bio(prices_test(i,:), q, costs, budget, elm_options, cnst_data, cnst_target) > 0) % constraint violation
            benefits_test(i) = 0;
        end
    end
    
    % [maxbenefits, maxbenefitsidx] = max(benefits_test);
    % fprintf(['      benefits:       £' sprintf('%s', num2sepstr(maxbenefits,  '%.0f')) '\n']);    
    % fprintf(['      budget surplus: £' sprintf('%s', num2sepstr(constraintfunc(prices_test(maxbenefitsidx,:)),'%.0f')) '\n']);    
    
    % Sort test_prices by benefits
    % ----------------------------
    [max_benefits_test, max_benefits_test_idx] = sort(benefits_test, 'descend');
    max_benefits_test_idx = max_benefits_test_idx(max_benefits_test > 0);
    
    % Return Nfeas prices
    % -------------------
    prices_feasible   = prices_test(max_benefits_test_idx(1:min(Nfeas,length(max_benefits_test_idx))),:);
    benefits_feasible = benefits_test(max_benefits_test_idx(1:min(Nfeas,length(max_benefits_test_idx))));
    
 end