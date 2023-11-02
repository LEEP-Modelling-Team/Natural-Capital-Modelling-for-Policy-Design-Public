% 1. Initialise
% -------------
clear
rng(23112010)

% Model
% -----
payment_mechanism = 'oc_pay';
unscaled_budget   = 1e9;
urban_pct_limit   = 0.5;
bio_constraint    = 0.15;
bio_as_prices     = false;
byparcel          = true;
carbon_price_string = 'non_trade_central';
drop_vars = {'habitat_non_use', 'biodiversity'};
budget_str  = [num2str(round(unscaled_budget/1e9)) 'bill'];
biocnst_str = [num2str(round(bio_constraint*100)) 'pct'];

% Markup
% ------
markup = 1.15;

% Paths to Data & Cplex Working Dir
% ---------------------------------
data_folder  = 'D:\myGitHub\defra-elms\Data\';
data_path = [data_folder, 'elm_data_', carbon_price_string, '.mat'];

% 2. Prepare data
% ---------------
data_year = 1;    
sample_size = 'no';  % all data
[b, c, q, hectares, budget, lu_data, cnst_data, cnst_target, elm_options, price_vars, new2kid] = ...
    load_data(sample_size, unscaled_budget, data_path, payment_mechanism, drop_vars, markup, urban_pct_limit, bio_constraint, bio_as_prices, byparcel, data_year);

% 3. Prepare data
% ---------------
% Create a table with column names
corr_tbl = table('Size', [0, 3], 'VariableNames', {'Option', 'Correlation', 'Value'}, 'VariableTypes', {'string', 'double', 'double'});

c_all = [];
b_all = [];

nproj = length(elm_options);
for i = 1:nproj
    
    % data for option type
    ci   = c(:,i);
    indi = ci>0;
    ci   = ci(indi);
    bi   = b(indi, i);

    % into c & b per hectare
    ci = table2array(ci./hectares(indi, i));
    bi = table2array(bi./hectares(indi, i));

    % correlation coefficient
    validind = ~isnan(ci) & ~isnan(bi);
    rho = corrcoef(ci(validind),bi(validind));
    val = mean(bi(validind));

    % store data
    newRow = {elm_options{i}, rho(1,2), val};
    corr_tbl = [corr_tbl; newRow];

end

