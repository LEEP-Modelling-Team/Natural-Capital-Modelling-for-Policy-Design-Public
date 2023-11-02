function f = myfun_ES(p, q, c, b, elm_options)

    num_farmers = length(c);
    num_options = length(elm_options);

    % Determine which option each farmer would prefer at these prices
    profit  = zeros(num_farmers, num_options+1);
    if isstruct(q) 
        for i = 1:num_options
            profit(:, i + 1)  = p * q.(elm_options{i})' - c(:, i)';
        end
    else
        for i = 1:num_options
            if ndims(q)>2
                profit(:, i + 1) = p * q(:, :, i)' - c(:, i)';
            else
                profit(:, i + 1) = p * q' - c(:, i)';
            end
        end
    end
    
	% Find which option gives each farmer maximum profit
	% They will choose 'do nothing' if all profits are negative
    [~, max_profit_col_idx] = max(profit, [], 2);
    
	% Calculate benefits for each farmer under this option uptake
    b = [zeros(num_farmers,1) b];
    b_chosen = b(sub2ind(size(b), (1:num_farmers)', max_profit_col_idx));
    
	% Calculate total benefits, return negative for minimisation
    f = -sum(b_chosen);
    
    % fprintf('%s\n', strjoin(num2sepstr(p,'% 3.5f'), ' '));
    % fprintf('%s\n', num2sepstr(f,'% 12.2f'));
    
end