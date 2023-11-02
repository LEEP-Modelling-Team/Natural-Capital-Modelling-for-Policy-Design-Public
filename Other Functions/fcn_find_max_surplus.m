
function max_surplus = fcn_find_max_surplus(q, c, prices_lb, prices_ub)

    num_cells   = size(q, 1);    
    num_prices  = size(q, 2);
    num_options = size(q, 3);

    surplus_min = zeros(num_cells, num_options);
    surplus_max = zeros(num_cells, num_options);
    
    for i = 1:num_options
        qi     = squeeze(q(:,:,i));
        pq_lb  = qi .* prices_lb';
        pq_ub  = qi .* prices_ub';
        pq_min = min(cat(3, pq_lb, pq_ub), [], 3);
        pq_max = max(cat(3, pq_lb, pq_ub), [], 3);

        rev_min = sum(pq_min,2);
        rev_max = sum(pq_max,2);

        surplus_min(:,i) = rev_min - c(:,i);
        surplus_max(:,i) = rev_max - c(:,i);
    end
    
    max_surplus = max(surplus_max, [], 2);
    max_surplus = surplus_max;
    
end