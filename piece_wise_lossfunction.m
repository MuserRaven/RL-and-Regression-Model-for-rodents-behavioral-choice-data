function [MAE] = piece_wise_lossfunction(alpha,r,beta,bias,gamma,mi_memo,choice_data,reward_data,choice_type_data, sound_cue_data)
% gamma_mi decay rate for mutual info
rng(10)
rat_choice = choice_data;
rat_reward = reward_data;  % reward matri
sound_cue = sound_cue_data;
choice_type = choice_type_data;

N = length(rat_choice); % timestep
q2=nan(N+1,1); %Expected Value for Left
q1=nan(N+1,1); %Exepected Value for Right
q2(1)=0; % Set initial Value to 0
q1(1)=0;

mutual_info1 = nan(N,1); % mutual information of choice 1
mutual_info2 = nan(N,1); % mutual information of choice 2
reward_dis2 = nan(N,1); % calculate info as entropy or simply estimated probability discrepancy of choice 2
choice1_percentage = length(rat_choice(rat_choice == 1 & choice_type == 1))/(length(rat_choice(rat_choice == 1 & choice_type == 1))+length(rat_choice(rat_choice == 2 & choice_type == 1)));
choice2_percentage = length(rat_choice(rat_choice == 2 & choice_type == 1))/(length(rat_choice(rat_choice == 1 & choice_type == 1))+length(rat_choice(rat_choice == 2 & choice_type == 1)));

for t=1:N

    % Compute mutual information 
    if isempty(find(sound_cue(1:t) == 1,1)) || isempty(find(sound_cue(1:t) == 2,1))
        reward_dis2(t) = 0; 
        mutual_info1(t) = 0; 
        mutual_info2(t) = 0; 

    else
    
            x1 = sound_cue(sound_cue(1:t) == 3 );
            y1 = rat_reward(sound_cue(1:t) == 3);
            mutual_info1(t) = mutualInfo(x1,y1); 
        
            x2 = sound_cue(sound_cue(1:t) == 1 | sound_cue(1:t) == 2 );
            y2 = rat_reward(sound_cue(1:t) == 1 | sound_cue(1:t) == 2);
            mutual_info2(t) = mutualInfo(x2,y2); 
    end
end

for t=1:N-1

  if choice_type(t) == 2  % remove force trials
        q1(t+1) =  q1(t) + alpha*(rat_reward(t) - q1(t)); % + r*mutual_info1(t)
        q2(t+1) =  (1-gamma)*(q2(t));

  elseif choice_type(t) == 3 
        q2(t+1) =  q2(t) + alpha*(rat_reward(t) - q2(t)); % + r*mutual_info2(t) 
        q1(t+1) =  (1-gamma)*(q1(t)) ;

  elseif rat_choice(t) == 1 && choice_type(t) == 1 % this means always no info 
    
        q1(t+1) =  q1(t) + alpha*(rat_reward(t) - q1(t)); % + r*mutual_info1(t);
        q2(t+1) =  (1-gamma)*(q2(t));
   
    
  elseif rat_choice(t) == 2 && choice_type(t) == 1
    
        q2(t+1) = q2(t) + alpha*(rat_reward(t) - q2(t)); %+ r*mutual_info2(t);
        q1(t+1) = (1-gamma)*(q1(t));
  end
end  

q1 = q1(1:N);  %(normalize(q1(1:N),'range',[0,1]));
q2 = q2(1:N); %(normalize(q2(1:N),'range',[0,1]));

%rat_choice(rat_choice == 1) = 0;rat_choice(rat_choice == 2) = 1;
%piece_wise_prob_choice2 = posterior_prob_choice(rat_choice);

sf = softmaxf(q2,q1,beta,bias);
sfni = softmaxf(q1,q2,beta,bias);

% Converting choice2_percentage into a binary value
choice2_binary = double(choice_data == 2);

% Calculate Cross-Entropy Loss



target_func = repmat(choice2_percentage,length(sf),1);
prob_info_corr = corr(mutual_info2,sf);

% Calculate MAE Loss
%MAE = mean(abs(sf(length(mutual_info2(mutual_info2 == 0)):end) - target_func(length(mutual_info2(mutual_info2 == 0)):end)));

% Calculate Cross-Entropy Loss
% Get indices where choice_type is free-choice
idx = (choice_type == 1);

% Subset your data based on these indices
choice2_binary_sub = choice2_binary(idx);
sf_sub = sf(idx);

% Calculate Mean Negative Log loss
epsilon = 1e-200;
MAE = -mean(choice2_binary_sub .* log(sf_sub + epsilon) + (1 - choice2_binary_sub) .* log(1 - sf_sub + epsilon));
