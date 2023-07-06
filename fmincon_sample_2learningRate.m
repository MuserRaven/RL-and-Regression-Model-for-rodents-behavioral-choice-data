%piece_wise_lossfunction(1,0.05,0.2,0.7,rat_choice,rat_reward,sound_cue,latency)
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
% baseline_data.CS(strcmp(baseline_data.CS, 'SPlus')) = {1};
% baseline_data.CS(strcmp(baseline_data.CS, 'SMinus')) = {2};
% baseline_data.CS(strcmp(baseline_data.CS, 'S3')) = {3};
id_list = unique(baseline_data.ID);
lambda = 0;% no regularization
test_condition = "ACC";
test_condition2 = "BLA";
test_condition3 = "OFC";
x = "F_Noinfo" ;
y = "F_Info" ;
corr_prob2_info2 = [];
suboptimal_info2 = []; % 1 = yes, 0 = 0
bias = 0;
count = 0;
nan_trial_time_rows = isnan(baseline_data.Trial_time);
nan_choice_rows = isnan(baseline_data.Choice);
nan_choice_rows_2 = isnan(baseline_data.Reinf);
nan_choice_rows_3 = isnan(baseline_data.RR);
% Find rows where either column has a NaN
rows_with_nan = nan_trial_time_rows | nan_choice_rows | nan_choice_rows_2 | nan_choice_rows_3;

% Use logical indexing to remove rows with NaN
baseline_data(rows_with_nan, :) = [];
% Create an empty table to store the parameters
% For 2 learing rate group
param_table_ACC_2lr = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'VariableNames', {'RatID','Sex','Expression','Group','Drug','Brain Area', 'Session', 'Alpha_r','Alpha_nr','R', 'Beta','bias','Gamma','MI_M(information decay)','Loss'});
param_table_control_2lr = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'VariableNames', {'RatID','Sex','Expression','Group','Drug','Brain Area', 'Session', 'Alpha_r','Alpha_nr','R', 'Beta','bias','Gamma','MI_M(information decay)','Loss'});
param_table_BLA_2lr = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'VariableNames', {'RatID','Sex','Expression','Group','Drug','Brain Area', 'Session', 'Alpha_r','Alpha_nr', 'R', 'Beta','bias','Gamma','MI_M(information decay)','Loss'});
param_table_OFC_2lr = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'VariableNames', {'RatID','Sex','Expression','Group','Drug','Brain Area', 'Session', 'Alpha_r','Alpha_nr', 'R', 'Beta','bias','Gamma','MI_M(information decay)','Loss'});


for patch = 1:4
    for id = 1:length(id_list)
        try
            rat_id = string(id_list{id});
            session_num = max(baseline_data.Session(baseline_data.ID == rat_id));

            concat_rat_choice = [];
            concat_rat_reward = [];
            concat_sound_cue = [];
            concat_choice_type = [];

            for session = 1:session_num
            if patch == 1
                trial_N_values = baseline_data.Trial_N(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");  
                [sorted_trial_N, sort_idx] = sort(trial_N_values);
                
                sex = baseline_data.Sex(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                drug = baseline_data.Drug(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                group = baseline_data.Group_s(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                expression = baseline_data.Group(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                brain_area = baseline_data.Brain(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_choice = baseline_data.Choice(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_reward = baseline_data.Reinf(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                choice_type = baseline_data.Trial_type(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                sound_cue = baseline_data.CS(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");

                rat_choice = rat_choice(sort_idx);
                rat_reward = rat_reward(sort_idx);
                choice_type = choice_type(sort_idx);
                sound_cue = sound_cue(sort_idx);
        

            elseif patch == 2
                trial_N_values = baseline_data.Trial_N(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                [sorted_trial_N, sort_idx] = sort(trial_N_values);
                
                sex = baseline_data.Sex(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                drug = baseline_data.Drug(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                group = baseline_data.Group_s(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                expression = baseline_data.Group(baseline_data.ID == rat_id & baseline_data.Group_s == "Control"& baseline_data.Session == session);
                brain_area = baseline_data.Brain(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                rat_choice = baseline_data.Choice(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                rat_reward = baseline_data.Reinf(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                choice_type = baseline_data.Trial_type(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);
                sound_cue = baseline_data.CS(baseline_data.ID == rat_id & baseline_data.Group_s == "Control" & baseline_data.Session == session);

                rat_choice = rat_choice(sort_idx);
                rat_reward = rat_reward(sort_idx);
                choice_type = choice_type(sort_idx);
                sound_cue = sound_cue(sort_idx);
        

             elseif patch == 3
                trial_N_values = baseline_data.Trial_N(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                [sorted_trial_N, sort_idx] = sort(trial_N_values);
                
                sex = baseline_data.Sex(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                drug = baseline_data.Drug(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                group = baseline_data.Group_s(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                expression = baseline_data.Group(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                brain_area = baseline_data.Brain(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_choice = baseline_data.Choice(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_reward = baseline_data.Reinf(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                choice_type = baseline_data.Trial_type(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                sound_cue = baseline_data.CS(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition2 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                
                rat_choice = rat_choice(sort_idx);
                rat_reward = rat_reward(sort_idx);
                choice_type = choice_type(sort_idx);
                sound_cue = sound_cue(sort_idx);
        

              elseif patch == 4
                trial_N_values = baseline_data.Trial_N(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                [sorted_trial_N, sort_idx] = sort(trial_N_values);
                
                sex = baseline_data.Sex(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                drug = baseline_data.Drug(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                group = baseline_data.Group_s(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                expression = baseline_data.Group(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                brain_area = baseline_data.Brain(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_choice = baseline_data.Choice(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                rat_reward = baseline_data.Reinf(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                choice_type = baseline_data.Trial_type(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");
                sound_cue = baseline_data.CS(baseline_data.ID == rat_id & baseline_data.Drug == "CNO"& baseline_data.Sex == "Female"& baseline_data.Brain == test_condition3 & baseline_data.Session == session & baseline_data.Group_s == "hM4Di");

                rat_choice = rat_choice(sort_idx);
                rat_reward = rat_reward(sort_idx);
                choice_type = choice_type(sort_idx);
                sound_cue = sound_cue(sort_idx);
        

            end
                % Concatenate rat_choice, rat_reward, and sound_cue across sessions
                concat_rat_choice = cat(1,concat_rat_choice, rat_choice);
                concat_rat_reward = cat(1,concat_rat_reward, rat_reward);
                concat_sound_cue = cat(1,concat_sound_cue,sound_cue);
                concat_choice_type = cat(1,concat_choice_type,choice_type);


            end


            rat_choice = concat_rat_choice;
            rat_reward = concat_rat_reward;
            sound_cue = concat_sound_cue;
            choice_type = concat_choice_type;

                        % Concatenating the vectors
            rat_data = [rat_choice, rat_reward, choice_type, sound_cue];
            
            % Converting matrix to table
            rat_data_table = array2table(rat_data, 'VariableNames', {'Choice', 'Reward', 'Trial_Type','Sound Cue'});
            
            % Removing rows with missing values
            %rat_data_table = rmmissing(rat_data_table);
    
           % Get the data
            rat_choice = table2array(rat_data_table(:,1));
            rat_reward = table2array(rat_data_table(:,2));
            choice_type = table2array(rat_data_table(:,3));
            sound_cue = table2array(rat_data_table(:,4));
            
            % Define objective function
            objective = @(x) piece_wise_lossfunction_2learningRate(x(1), x(2), 0, x(3),x(4), x(5),inf, rat_choice, rat_reward, choice_type, sound_cue);
            
            
            % Lower and upper bounds for alpha, r, and beta
            % alpha_r,alpha_nr,r,beta,bias,gamma,gamma_mi,
            lb = [eps,eps, eps, -1000,eps];  % change to the lower bounds suitable for your problem
            ub = [1, 1, 100, 1000,1];  % change to the upper bounds suitable for your problem
            x0 = (lb + ub)/2;  % Midpoint of bounds
            
            % Create a problem structure
            problem = createOptimProblem('fmincon','objective',objective,'x0',x0,'lb',lb,'ub',ub);

           
            % Define the multi-start object
            ms = MultiStart;
            
            % Number of multiple starting points
            nStarts = 500;  % change this to the desired number of starting points
            
            % Perform the optimization
            [x,fval] = run(ms,problem,nStarts);

            % After optimization, compute the loss and add results to the correct table
            loss =  piece_wise_lossfunction_2learningRate(x(1), x(2), 0, x(3),x(4), x(5),inf,rat_choice,rat_reward, choice_type, sound_cue);

            if patch == 1
                param_table_ACC_2lr = [param_table_ACC_2lr; {rat_id,sex, expression,group,drug, brain_area, session, x(1), x(2), 0, x(3),x(4), x(5),inf,loss}];
            elseif patch == 2
                param_table_control_2lr = [param_table_control_2lr; {rat_id,sex, expression,group,drug, brain_area, session, x(1), x(2), 0, x(3),x(4), x(5),inf, loss}];
            elseif patch == 3
                param_table_BLA_2lr = [param_table_BLA_2lr; {rat_id,sex, expression,group,drug, brain_area, session, x(1), x(2), 0, x(3),x(4), x(5),inf,loss}];
            elseif patch == 4
                param_table_OFC_2lr = [param_table_OFC_2lr; {rat_id,sex, expression,group,drug, brain_area, session, x(1), x(2), 0, x(3),x(4), x(5),inf,loss}];
            end

            count = count + 1;
        catch ME % Catch the error

            fprintf('Error in id %d: %s\n', id, ME.message); % Display error message
            continue; % Continue with the next iteration
        end
    end
end
%% Log Loss compute

% Extract "Loss" column from each table
loss_ACC = param_table_ACC_2lr.Loss;
loss_BLA = param_table_BLA_2lr.Loss;
loss_control = param_table_control_2lr.Loss;
loss_OFC = param_table_OFC_2lr.Loss;

% Concatenate the data from all tables
loss_all = [loss_ACC; loss_BLA; loss_control; loss_OFC];

% Create a new figure
figure;

% Plot the histogram with more bins
num_bins = 15; % Change this value to modify the number of bins
histogram(loss_all, num_bins, 'FaceColor', [0.7 0.7 0.7],'Normalization','pdf');hold on
xline(mean(loss_all),'LineWidth',3,'LineStyle','--','Color','r');hold on
text(mean(loss_all)-0.2,0.8,{'Mean = ', mean(loss_all)})
xline(0.693,'LineWidth',3,'LineStyle','--','Color','b');
legend('','mean loss','baseline model loss')
length(find(loss_all<0.693))/length(loss_all)
% Set the title and labels
title('Mean NLL Loss across all conditions', 'FontSize', 14);
xlabel('Loss', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);

% Set the background color of the figure to white
set(gcf, 'color', 'w');

% Set axes properties
set(gca, 'LineWidth', 2, 'FontSize', 12, 'Box', 'off');

%% Correlation Plot
% Define the columns of interest (MATLAB uses 1-indexing)
cols_of_interest = 8:12;

% Extract the columns of interest from each table
ACC_data = param_table_ACC_2lr(:, cols_of_interest);
BLA_data = param_table_BLA_2lr(:, cols_of_interest);
OFC_data = param_table_OFC_2lr(:, cols_of_interest);
Control_data = param_table_control_2lr(:, cols_of_interest);

% Compute the correlation matrices
ACC_corr = corr(table2array(ACC_data), 'Type', 'Pearson');
BLA_corr = corr(table2array(BLA_data), 'Type', 'Pearson');
OFC_corr = corr(table2array(OFC_data), 'Type', 'Pearson');
Control_corr = corr(table2array(Control_data), 'Type', 'Pearson');

% Create the correlation heatmaps for each condition
f

subplot(2,2,1);
heatmap(ACC_data.Properties.VariableNames, ACC_data.Properties.VariableNames, ACC_corr, 'Colormap', jet, 'MissingDataColor', 'white');
title('ACC Pairwise Correlation');

subplot(2,2,2);
heatmap(BLA_data.Properties.VariableNames, BLA_data.Properties.VariableNames, BLA_corr, 'Colormap', jet, 'MissingDataColor', 'white');
title('BLA Pairwise Correlation');

subplot(2,2,3);
heatmap(OFC_data.Properties.VariableNames, OFC_data.Properties.VariableNames, OFC_corr, 'Colormap', jet, 'MissingDataColor', 'white');
title('OFC Pairwise Correlation');

subplot(2,2,4);
heatmap(Control_data.Properties.VariableNames, Control_data.Properties.VariableNames, Control_corr, 'Colormap', jet, 'MissingDataColor', 'white');
title('Control Pairwise Correlation');

% Create the colorbar
cb = colorbar('Position', [0.93 0.11 0.02 0.815]);
ylabel(cb, 'Pearson Correlation Coefficient');
%% BOx plot

tables = {param_table_ACC_2lr, param_table_BLA_2lr, param_table_OFC_2lr, param_table_control_2lr};
labels = {'ACC', 'BLA', 'OFC', 'Control'};
% Number of bootstrap samples
n_bootstrap = 1000;
% Define the columns of interest (MATLAB uses 1-indexing)
cols_of_interest = [8,9,11,12,13];

% Create a new figure
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]); % you can adjust this as per your requirements

% Color scheme for boxplots
colors = {[1, 0, 0], [0, 0, 1], [0, 0.75, 0.75], [0, 1, 0], [0, 1, 0]};  % Adjusted color for 'Control' condition to green

for j = 1:numel(cols_of_interest)
    
    % Using subplot for flexibility
    subplot(1, numel(cols_of_interest), j);
    hold on;

    bs_samples = cell(numel(tables), 1);
    for i = 1:numel(tables)
        data = tables{i}(:, cols_of_interest(j)).Variables;
        bs_mean = [];
        bs_cv = [];
        bs_25 = [];
        bs_75 = [];
        bs_2_5 = []; % 2.5 percentile for 95% CI
        bs_97_5 = []; % 97.5 percentile for 95% CI
        for h = 1:n_bootstrap
            sample = datasample(data, length(data));
            bs_mean = cat(1,bs_mean,median(sample));
            bs_cv = cat(1,bs_cv,std(sample)/mean(sample));  % calculate coefficient of variation
            bs_25 = cat(1,bs_25,prctile(sample, 12.5));
            bs_75 = cat(1,bs_75,prctile(sample, 87.5));
            bs_2_5 = cat(1,bs_2_5,prctile(sample, 2.5));
            bs_97_5 = cat(1,bs_97_5,prctile(sample, 97.5));
        end

        bs_mean = mean(bs_mean); 
        bs_cv = mean(bs_cv);
        bs_25 = mean(bs_25);
        bs_75 = mean(bs_75);
        bs_2_5 = mean(bs_2_5);
        bs_97_5 = mean(bs_97_5);
        
         % Printing out the bootstrap statistics
        fprintf('Condition: %s, Parameter: %s\n', labels{i}, tables{i}.Properties.VariableNames{cols_of_interest(j)});
        fprintf('Bootstrap Mean: %.2f\n', bs_mean);
        fprintf('Bootstrap Coeff. of Variation: %.2f\n', bs_cv);
        fprintf('Bootstrap 25th percentile: %.2f\n', bs_25);
        fprintf('Bootstrap 75th percentile: %.2f\n', bs_75);
        fprintf('Bootstrap 2.5th percentile (95%% CI Lower): %.2f\n', bs_2_5);
        fprintf('Bootstrap 97.5th percentile (95%% CI Upper): %.2f\n\n', bs_97_5);

        % Error bar for 95% CI
        line([i i], [bs_2_5 bs_97_5], 'Color', 'k', 'LineWidth', 2);  % vertical line
        line([i-0.15 i+0.15], [bs_2_5 bs_2_5], 'Color', 'k', 'LineWidth', 2);  % bottom horizontal line
        line([i-0.15 i+0.15], [bs_97_5 bs_97_5], 'Color', 'k', 'LineWidth', 2);  % top horizontal line

        % Box for 25th-75th percentiles with semi-transparent color
        rectangle('Position',[i-0.25 bs_25 0.5 bs_75-bs_25], 'FaceColor', [colors{i}, 1], 'EdgeColor', 'k'); 
        % Line for mean
        line([i-0.25 i+0.25], [bs_mean bs_mean],'Color', 'k', 'LineWidth', 2);  

    end
    
    % Set the title based on the index
    if j == 1
        title('Alpha r', 'FontSize', 14, 'FontWeight', 'bold');
    elseif j == 2
        title('Alpha nr', 'FontSize', 14, 'FontWeight', 'bold');
    else
        title(tables{1}.Properties.VariableNames{cols_of_interest(j)}, 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    if j == 1
        ylabel('Value', 'FontSize', 12, 'FontWeight', 'bold');  % Only for the first subplot
    end
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'FontWeight', 'bold');
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2, 'FontSize', 12, 'Box', 'off');
       % Add 'Condition' label at the bottom of each subplot
%     xlabel('Condition', 'FontSize', 12, 'FontWeight', 'bold');
  
    hold off;
end
%% No bootsrap

tables = {param_table_ACC_2lr, param_table_BLA_2lr, param_table_OFC_2lr, param_table_control_2lr};
labels = {'ACC', 'BLA', 'OFC', 'Control'};

cols_of_interest = [8,9,11,12,13];

% Define colors for boxplots
colors = {'r', 'b', 'g', 'k'};

% Define positions for subplots
positions = {[0.1 0.6 0.4 0.3], [0.6 0.6 0.4 0.3], [0.1 0.1 0.4 0.3], [0.6 0.1 0.4 0.3]};

figure
bs_means = cell(numel(tables), numel(cols_of_interest)); 
bs_cvs = cell(numel(tables), numel(cols_of_interest));  

for j = 1:numel(cols_of_interest)
    subplot('Position', positions{j});
    hold on;
    bs_samples = cell(numel(tables), 1);
    for i = 1:numel(tables)
        data = tables{i}(:, cols_of_interest(j)).Variables;
        bs_samples{i} = datasample(data, n_bootstrap);
        bs_means{i,j} = mean(bs_samples{i});
        bs_cvs{i,j} = std(bs_samples{i}) / bs_means{i,j};
        fprintf('Bootstrap mean and cv for Condition %s, Feature %s: %.2f, %.2f\n', labels{i}, tables{i}.Properties.VariableNames{cols_of_interest(j)}, bs_means{i,j}, bs_cvs{i,j});
        h = boxplot(data, 'Colors', colors{i}, 'Widths', 0.5, 'Positions', i, 'Symbol', '');
        set(h, 'linewidth', 2);
    end
    title(tables{1}.Properties.VariableNames{cols_of_interest(j)}, 'FontSize', 14);
    ylabel('Value', 'FontSize', 12);
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2, 'FontSize', 12, 'Box', 'off');
    hold off;
end