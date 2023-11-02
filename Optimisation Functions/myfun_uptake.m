function uptake = myfun_uptake(p, q, c, elm_option)

    num_farmers = length(c);
    num_options = length(elm_option);

    % Determine which option each farmer would prefer at these prices
    profit  = zeros(num_farmers, num_options+1);
    
    % Evaluate profit & cost for each option 
    if isstruct(q) 
        for i = 1:num_options
            profit(:, i + 1)  = p * q.(elm_option{i})' - c(:, i)';
        end     
    else
        if ndims(q) == 2
             profit(:, 2:num_options+1) = p .* q - c;
        else
            for i = 1:num_options
                profit(:, i + 1) = p * q(:, :, i)' - c(:, i)';
            end
        end
    end
    
	% Find which option gives each farmer maximum profit
	% They will choose 'do nothing' if all profits are negative
    [~, max_profit_col_idx] = max(profit, [], 2);

    % uptake
    uptake = full(sparse(1:num_farmers, max_profit_col_idx, ones(num_farmers,1), num_farmers, num_options+1)); 
    
    % remove do-nothing option
    uptake = uptake(:,2:end);
    
end