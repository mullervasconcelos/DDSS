clear
clc

global n xhat data theta

load testing_data.mat

data1 = [T_out T1 T2 T3 T4 T5 T6 T7 T8 T9 RH_out RH_1 RH_2 RH_3 RH_4 RH_5 RH_6 RH_7 RH_8 RH_9];


mean_data1 = mean(data1)

% Remove a few missing data
data1(1,:) = [];

mean_data1 = mean(data1)

load training_data.mat

data2 = [T_out T1 T2 T3 T4 T5 T6 T7 T8 T9 RH_out RH_1 RH_2 RH_3 RH_4 RH_5 RH_6 RH_7 RH_8 RH_9];

mean_data2 = mean(data2)

% Remove a few missing data
data2(1,:) = [];

mean_data2 = mean(data2)

data = [data1; data2];



[N,n] = size(data)