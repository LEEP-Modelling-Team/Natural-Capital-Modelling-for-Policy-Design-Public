function spend = myfun_spend(p, q, c, elm_options)

    num_farmers = length(c);
    num_options = length(elm_options);

    % Determine which option each farmer would prefer at these prices
    profit = zeros(num_farmers, num_options+1);
    spend  = zeros(num_farmers, num_options);
    % Evaluate profit & cost for each option 
    if isstruct(q) 
        for i = 1:num_options
            payment       = p * q.(elm_options{i})';
            profit(:,i+1) = payment - c(:, i)';
            spend(:,i)    = payment;               
        end     
    else
        if ndims(q) == 2
            spend = p .* q;
            profit(:, 2:num_options+1) = spend - c;
        else
            for i = 1:num_options
                payment       = p * q(:, :, i)';
                profit(:,i+1) = payment - c(:, i)';
                spend(:,i)    = payment;                   
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
    
    spend = sum(uptake.*spend,'all');
    
end