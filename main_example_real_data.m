clear
clc

global n xhat data theta

load training_data.mat

data_train = [T10 RH_10 T11 RH_11 T12 RH_12 T13 RH_13 T14 RH_14 T16 RH_16 T17 RH_17 T18 RH_18];

% Remove a few missing data
data_train(1,:) = [];

mean_data1 = mean(data_train)

[N1,n1] = size(data_train)

load testing_data.mat

data_test = [T1 RH_1 T2 RH_2 T3 RH_3 T4 RH_4 T5 RH_5 T7 RH_7 T8 RH_8 T9 RH_9];
    
% Remove a few missing data
data_test(1,:) = [];

mean_data2 = mean(data_test)

[N2,n2] = size(data_test)

data = [data_train; data_test];

%hist(data(:,1),100)

[N,n] = size(data);

p = randperm(N);

data_mixed = data(p,:);

s = floor(N/21);

train_batch = data_mixed(1:s,:);

test_batches = data_mixed(s+1:s+20*s,:);

size(test_batches);

save('train_batch.mat','train_batch')

save('test_batches.mat','test_batches')
