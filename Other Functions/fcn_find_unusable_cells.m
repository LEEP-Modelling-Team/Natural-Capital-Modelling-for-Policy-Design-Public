function unusable_cells = fcn_find_unusable_cells(q, c, budget, elm_options, prices_max, prices)
    
    % Constants
    % ---------
    num_farmers = size(c,1);
    num_options = size(c,2);
    num_prices  = length(prices_max);

    % Choose some test prices
    % -----------------------
    % o prices_max are prices if spent whole budget on just buying the one
    %   service. Put on diagonal so test price vector is that price and 
    %   zero for all other services.
    % o prices are a set of additional prices where want to make sure that
    %   all farms that uptake are those prices are in the sample.

    prices_test = diag(prices_max);
    if exist('prices', 'var')
       prices_test = [prices_test; prices'];
    end
    prices_test(isinf(prices_test)) = 0;

    % Create combinations of max prices using latin hyercube
    % ------------------------------------------------------
    N = 1000;
    prices_test = [prices_test; lhsdesign(N, num_prices).*prices_max'];
    
    % Exclude if Farm not Used at any of test prices
    % ----------------------------------------------
    uptake = zeros(num_farmers, num_options);
    for ii = 1:length(prices_test)
        uptake = uptake + myfun_uptake(prices_test(ii,:), q, c, elm_options);
        % fprintf('%0.0f. Num chosen: %0.0f \n', ii, sum(uptake>0, 'all'));
    end
    % fprintf('%0.0f. Num chosen: %0.0f \n', ii, sum(uptake>0, 'all'));
    unusable_cells = sum(uptake, 2) == 0; 

end