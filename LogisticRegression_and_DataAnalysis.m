%% Run Logistic Regression



% Define the objective function
rng(10)
% Preprocessing the file:
% % For the most recent data:
steady_state_data = readtable('1_Steadystate_R.xlsx');
% 
% % Convert cell array to categorical array
steady_state_data.CS = categorical(steady_state_data.CS);
steady_state_data.Trial_type = categorical(steady_state_data.Trial_type);
% 
% % Replace values in 'CS' column
CS_numeric = zeros(size(steady_state_data.CS));
CS_numeric(steady_state_data.CS == 'SPlus') = 1;
CS_numeric(steady_state_data.CS == 'SMinus') = 2;
CS_numeric(steady_state_data.CS == 'S3') = 3;
steady_state_data.CS = CS_numeric;
% 
% % Replace values in 'Trial_type' column
Trial_type_numeric = zeros(size(steady_state_data.Trial_type));
Trial_type_numeric(steady_state_data.Trial_type == 'Choice') = 1;
Trial_type_numeric(steady_state_data.Trial_type == 'F_Noinfo') = 2;
Trial_type_numeric(steady_state_data.Trial_type == 'F_Info') = 3;
steady_state_data.Trial_type = Trial_type_numeric;



baseline_data = steady_state_data;
id_list = unique(baseline_data.ID);



% Initialize empty table to store p-values
%pValueTable_control = table();
%pValueTable_test_ACC = table();
pValueTable_test_BLA = table();
%pValueTable_test_OFC = table();
count = 0; 
R2_control = 0;
oddsRatiosTable_control = table();
ciTable_control = table();
ciTable_control80 = table();
oddsRatiosTable_ACC = table();
ciTable_ACC = table();
ciTable_ACC80 = table();
oddsRatiosTable_BLA = table();
ciTable_BLA = table();
ciTable_BLA80 = table();
oddsRatiosTable_OFC = table();
ciTable_OFC = table();
ciTable_OFC80 = table();

id_list = unique(baseline_data.ID);
% Loop over each subject
for patch = 1:4
% for i = 1:length(id_list)
% 
%     rat_id = string(id_list{i});
if patch == 1
    rat_choice = baseline_data.Choice( baseline_data.Drug == "CNO" & baseline_data.Sex == "Female"& baseline_data.Brain == "ACC" & baseline_data.Group_s == "hM4Di");
    concat_choice_type = baseline_data.Trial_type(baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "ACC" & baseline_data.Group_s == "hM4Di");
    rat_reward = baseline_data.Reinf(baseline_data.Drug== "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "ACC" & baseline_data.Group_s == "hM4Di");
    rat_latency = baseline_data.Lat_Choice( baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "ACC" & baseline_data.Group_s == "hM4Di");
elseif patch == 2  
    rat_choice = baseline_data.Choice( baseline_data.Group_s == "Control" );
    concat_choice_type = baseline_data.Trial_type(baseline_data.Group_s == "Control" );
    rat_reward = baseline_data.Reinf( baseline_data.Group_s == "Control" );
    rat_latency = baseline_data.Lat_Choice(baseline_data.Group_s == "Control" );
elseif patch == 3
    rat_choice = baseline_data.Choice(baseline_data.Drug == "CNO" & baseline_data.Sex == "Female"& baseline_data.Brain == "BLA" & baseline_data.Group_s == "hM4Di");
    concat_choice_type = baseline_data.Trial_type(baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "BLA" & baseline_data.Group_s == "hM4Di");
    rat_reward = baseline_data.Reinf(baseline_data.Drug== "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "BLA" & baseline_data.Group_s == "hM4Di");
    rat_latency = baseline_data.Lat_Choice(baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "BLA" & baseline_data.Group_s == "hM4Di");
elseif patch == 4
    rat_choice = baseline_data.Choice(baseline_data.Drug == "CNO" & baseline_data.Sex == "Female"& baseline_data.Brain == "OFC" & baseline_data.Group_s == "hM4Di");
    concat_choice_type = baseline_data.Trial_type(baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "OFC" & baseline_data.Group_s == "hM4Di");
    rat_reward = baseline_data.Reinf(baseline_data.Drug== "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "OFC" & baseline_data.Group_s == "hM4Di");
    rat_latency = baseline_data.Lat_Choice(baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == "OFC" & baseline_data.Group_s == "hM4Di");
end

    if isempty(rat_choice)
        continue
    end


    % Concatenating the vectors
    rat_data = [rat_choice, rat_reward, concat_choice_type, rat_latency];
    
    % Converting matrix to table
    rat_data_table = array2table(rat_data, 'VariableNames', {'Choice', 'Reward', 'Trial_Type','Latency'});

    rat_choice = table2array(rat_data_table(:,1));
    rat_reward = table2array(rat_data_table(:,2));
    choice_type = table2array(rat_data_table(:,3));
    rat_latency = table2array(rat_data_table(:,4));


    % Adjusting rat_choice to have the values 0 and 1
    rat_choice(rat_choice == 1) = 0;
    rat_choice(rat_choice == 2) = 1;
    
    % Find indices where concat_choice_type == 1
    indices = find(concat_choice_type == 1);
    
    % Shift indices 1 step before
    indices_shifted = indices - 1;
    filtered_rat_choice = rat_choice(indices);
    
    % Ensure no index is less than 1
    indices_shifted(indices_shifted < 1) = 1;
    
    % Filter rat_choice and rat_reward
    previous_choice = rat_choice(indices_shifted);
    previous_reward = rat_reward(indices_shifted);
    filtered_latency = rat_latency(indices);
    %length(filtered_rat_choice)

    % Create table for the logistic regression
    rat_data = table(filtered_rat_choice(2:end), previous_choice(2:end), previous_reward(2:end), filtered_latency(2:end), ...
        'VariableNames', {'Choice', 'PrevChoice', 'PrevReward', 'Latency'});
    count = count +1;
    % Create a new row
    % Define options for maximum iteration
    options = statset('MaxIter', 1000000000000, 'TolFun', 1e-6, 'TolX', 1e-6);
    mdl = fitglm(rat_data, 'Choice ~ PrevChoice + PrevReward + Latency', 'Distribution', 'binomial', 'Options', options);

    % Calculation of the odds ratio
    odds_ratios = exp(mdl.Coefficients.Estimate);

    % Create new rows for odds ratio and CI tables
    newRow_oddsRatios = array2table(odds_ratios', 'VariableNames', {'Intercept', 'PrevChoice', 'PrevReward', 'Latency'});
         
    ci_lower = exp(mdl.Coefficients.Estimate - 1.96 * mdl.Coefficients.SE);
    ci_upper = exp(mdl.Coefficients.Estimate + 1.96 * mdl.Coefficients.SE);
    % Combine ci_lower and ci_upper into one matrix
    combined_ci = [ci_lower, ci_upper];
    
    % Reshape the combined_ci into a single row
    reshaped_ci = reshape(combined_ci', 1, []);
    
    % Create a table
    newRow_ci = array2table(reshaped_ci, 'VariableNames', ...
                            {'Intercept_lower', 'Intercept_upper', 'PrevChoice_lower', 'PrevChoice_upper', 'PrevReward_lower', 'PrevReward_upper', 'Latency_lower', 'Latency_upper'});
    

    ci_lower80 = exp(mdl.Coefficients.Estimate - 1.28 * mdl.Coefficients.SE);
    ci_upper80 = exp(mdl.Coefficients.Estimate + 1.28 * mdl.Coefficients.SE);
    % Combine ci_lower and ci_upper into one matrix
    combined_ci80 = [ci_lower80, ci_upper80];
    
    % Reshape the combined_ci into a single row
    reshaped_ci80 = reshape(combined_ci80', 1, []);
    
    % Create a table
    newRow_ci80 = array2table(reshaped_ci80, 'VariableNames', ...
                            {'Intercept_lower', 'Intercept_upper', 'PrevChoice_lower', 'PrevChoice_upper', 'PrevReward_lower', 'PrevReward_upper', 'Latency_lower', 'Latency_upper'});
    


    % Append the new rows to the odds ratio and CI tables
    if patch == 1
        oddsRatiosTable_ACC = [oddsRatiosTable_ACC; newRow_oddsRatios];
        ciTable_ACC = [ciTable_ACC; newRow_ci];
        ciTable_ACC80 = [ciTable_ACC80; newRow_ci80];
    elseif patch ==2
        oddsRatiosTable_control = [oddsRatiosTable_control; newRow_oddsRatios];
        ciTable_control = [ciTable_control; newRow_ci];
        ciTable_control80 = [ciTable_control80; newRow_ci80];
    elseif patch == 3
        oddsRatiosTable_BLA = [oddsRatiosTable_BLA; newRow_oddsRatios];
        ciTable_BLA = [ciTable_BLA; newRow_ci];
        ciTable_BLA80 = [ciTable_BLA80; newRow_ci80];
    elseif patch == 4
        oddsRatiosTable_OFC = [oddsRatiosTable_OFC; newRow_oddsRatios];
        ciTable_OFC = [ciTable_OFC; newRow_ci];
        ciTable_OFC80 = [ciTable_OFC80; newRow_ci80];
    end

%     if mdl.Rsquared.Ordinary ~= -inf
%         R2_control = mdl.Rsquared.Ordinary;
%     else
%         R2_control = R2_control + 0;
%         count = count - 1;
%     end
    mdl.Rsquared.Adjusted
%end
end
%% Logistic Regression Data Analysis
% Create figure and axis handles
figure;

% Define labels and colors for conditions
conditionLabels = {'ACC', 'BLA', 'OFC', 'Control'};
colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1]}; % Define colors for conditions (red, green, blue, cyan)

% Store tables in cell arrays for easier loop access
oddsRatioTables = {oddsRatiosTable_ACC, oddsRatiosTable_BLA, oddsRatiosTable_OFC, oddsRatiosTable_control};
ciTables80 = {ciTable_ACC80, ciTable_BLA80, ciTable_OFC80, ciTable_control80};
ciTables95 = {ciTable_ACC, ciTable_BLA, ciTable_OFC,ciTable_control};

% Loop over features
for featureIdx = 2:4  % columns 2:4 correspond to 'PrevChoice', 'PrevReward', 'Latency'
    subplot(1, 3, featureIdx - 1);
    hold on;
    
    % Loop over conditions
    for condIdx = 1:4
        % Extract corresponding odds ratios and CIs
        oddsRatio = oddsRatioTables{condIdx}{:, featureIdx};
        ci80 = ciTables80{condIdx}{:, 2*(featureIdx - 1) + (1:2)};  % extract lower and upper CI (80%)
        ci95 = ciTables95{condIdx}{:, 2*(featureIdx - 1) + (1:2)};  % extract lower and upper CI (95%)

        % Plot
        % Note: you might want to adjust positions (e.g., condIdx) and/or add scaling factors
        % to better separate overlapping elements in the plot
        errorbar(condIdx, oddsRatio, oddsRatio - ci95(1), ci95(2) - oddsRatio, 'k', 'LineWidth', 2);  % CI 95%
        rectangle('Position', [condIdx - 0.25, ci80(1), 0.5, ci80(2) - ci80(1)], 'FaceColor', colors{condIdx});  % CI 80%
        line([condIdx - 0.25, condIdx + 0.25], [oddsRatio, oddsRatio], 'Color', 'k', 'LineWidth', 2);  % Odds Ratio
    end
    
    % Add some labels, title, legend, etc.
    set(gca, 'XTick', 1:4, 'XTickLabel', conditionLabels, 'LineWidth', 2);
    xlabel('Condition');
    ylabel('Odds Ratio');
    title([ oddsRatioTables{1}.Properties.VariableNames{featureIdx}]);
    hold off;
end
