function [c, ceq] = mycon_budget(p, q, cst, budget, elm_options)

    num_farmers = length(cst);
    num_options = length(elm_options);

    % Determine which option each farmer would prefer at these prices
    profit = zeros(num_farmers, num_options + 1);
    spend  = zeros(num_farmers, num_options + 1);
	
    % Evaluate profit and spend for each option and farmer
    if isstruct(q) 
        for i = 1:num_options
            profit(:, i + 1) = p * q.(elm_options{i})' - cst(:, i)';
            spend(:, i + 1)  = p * q.(elm_options{i})';
        end
    elseif ndims(q) > 2
        for i = 1:num_options
            profit(:, i + 1) = p * q(:, :, i)' - cst(:, i)';
            spend(:, i + 1)  = p * q(:, :, i)';
        end
    else
        profit(:, 2:end) = p .* q - cst;
        spend(:, 2:end)  = p .* q;        
    end
    
	% Find which option gives each farmer maximum profit
	% They will choose 'do nothing' if all profits are negative
    [~, max_profit_col_idx] = max(profit, [], 2);
    
	% Calculate spend for each farmer under this option uptake
    spend_final = spend(sub2ind(size(spend), (1:num_farmers)', max_profit_col_idx));
    
	% Calculate margin (total spend - budget)
    c   = sum(spend_final) - budget;
    ceq = [];
    
end