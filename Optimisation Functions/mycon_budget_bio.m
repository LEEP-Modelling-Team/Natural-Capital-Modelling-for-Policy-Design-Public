function [c, ceq] = mycon_budget_bio(p, q, costs, budget, elm_options, cnst_data, cnst_target)

    % Constants
    % --------
    N  = size(costs,1);
    Nq = length(elm_options);

    % Calculate Profit & Spend
    % ------------------------ 
    % Determine which option each farmer would prefer at these prices
    profit = zeros(N, Nq + 1);
    spend  = zeros(N, Nq + 1);	
    % Evaluate profit and spend for each option and farmer
    if isstruct(q) 
        for i = 1:Nq
            profit(:, i + 1) = p * q.(elm_options{i})' - costs(:, i)';
            spend(:, i + 1)  = p * q.(elm_options{i})';
        end
    else
        for i = 1:Nq
            if ndims(q)>2
                profit(:, i + 1) = p * q(:, :, i)' - costs(:, i)';
                spend(:, i + 1)  = p * q(:, :, i)';
            else
                profit(:, i + 1) = p * q' - costs(:, i)';
                spend(:, i + 1)  = p * q';
            end                        
        end
    end
    
    % Calculate Uptake Matrix
    % -----------------------        
	% Find which option gives each farmer maximum profit
	% They will choose 'do nothing' if all profits are negative
    [~, max_profit_col_idx] = max(profit, [], 2);
    % uptake
    uptake = full(sparse(1:N, max_profit_col_idx, ones(N,1), N, Nq+1));     
    % remove do-nothing option
    uptake = uptake(:,2:end);
        
    % Check Budget Constraint
    % -----------------------    
	% Calculate spend for each farmer under this option uptake
    spend_final = sum(spend(:,2:end) .* uptake, 2);
    
    % Budget constraint: spend - budget < 0
    c   = sum(spend_final) - budget;
    ceq = [];
    
    % Check Biodiversity Constraint
    % -----------------------------
    if any(cnst_target)

        num_spgrp = length(cnst_target);
        
        spgrp_chg = zeros(num_spgrp,1);
        for k = 1:num_spgrp
            spgrp_chg(k) = sum(uptake.*squeeze(cnst_data(k,:,:))', 'all');        
        end
        
        % Biodiversity target constraints: cnst_target - spgrp_chg <0
        c = [c; cnst_target - spgrp_chg];
        
    end
        
end