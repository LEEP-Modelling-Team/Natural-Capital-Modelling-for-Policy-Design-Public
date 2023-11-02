function prices_max = fcn_find_max_prices(q, c, budget, lu_data, elm_options)

    % Constants    
    % ---------
    num_farmers = size(q,1);
    num_prices  = size(q,2); 
    num_options = size(q,3);
    prices_max  = zeros(num_prices, 1);
    
    % cost to infinite for options not for parcel land use
    c(lu_data==0) = inf;
    
    % Constants    
    % ---------   
    
    for j = 1:num_prices        
        qj = squeeze(q(:,j,:));    
        price_pts = c./qj;
        price_pts(isinf(price_pts)) = NaN;
        price_pts = sort(price_pts(:), 'ascend', 'MissingPlacement', 'last');
               
        % Binary criterion function
        inbudget = @(p) myfun_spend(p,qj,c,elm_options) < budget;

        % Initialize the high and low bounds of the search
        low = 1;
        high = length(price_pts);
        
        % Initialize the highest value that meets the criterion
        highest = NaN;

        % Perform a binary search until the high and low bounds converge
        while low <= high
            % Calculate the middle index and value of the search range
            mid   = floor((low+high)/2);
            price = price_pts(mid);

            % Check if the value meets the criterion
            if inbudget(price)
                % If it does, update the highest value that meets the criterion
                highest = price;

                % Search the upper half of the range
                low = mid+1;
            else
                % If it doesn't, search the lower half of the range                    
                if isnan(price)||price > 0
                    high = mid-1;
                else
                    % unless price was negative and then got to go back up
                    % again
                    low = mid+1;
                end                      
            end
        end  
        
        prices_max(j) = highest;
        
    end
    
end    
    