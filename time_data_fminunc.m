
readtable('Time_data_fminunc_unicast.xlsx')

time_data = [ans.T1 ans.T2 ans.T3 ans.T4 ans.T5 ans.T6 ans.T7 ans.T8 ans.T9 ans.T10]';

n_var = [2 5 10 20 50 100]

errorbar(n_var, mean(time_data), std(time_data));

hold

readtable('Time_data_ccp_unicast.xlsx')

time_data = [ans.T1 ans.T2 ans.T3 ans.T4 ans.T5 ans.T6 ans.T7 ans.T8 ans.T9 ans.T10]';

n_var = [2 5 10 20 50 100]

errorbar(n_var, mean(time_data), std(time_data));