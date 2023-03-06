function good = Find_Period(t, Dim5, diff_base, var_base, mean_base) % history function for t < t0
global TAU Var Diff window_Time window_TAU
good=0;
delt_t = t(2) - t(1);
TAU = [];
Var = [];
Diff = [];
window_TAU = [];
window_Time = [];
sum_tau = [Dim5(end-1, 1)];       % sum_tau[i] = tau[0]   + ... + tau[i]
sum_tau2 = [Dim5(end-1, 1) ^ 2]; % sum_tau[i] = tau[0]^2 + ... + tau[i]^2
min_diff = 100000;

for i=2:length(t)
    t_now = t(i);
    now_tau = Dim5(end-1, i);
    sum_tau = [sum_tau, sum_tau(i - 1) + now_tau];
    sum_tau2 = [sum_tau2, sum_tau2(i - 1) + now_tau ^ 2];
    interval_len = fix(now_tau / delt_t);
    if interval_len + 5 > i
        continue
    end
    mean_tau = (sum_tau(i) - sum_tau(i - interval_len - 1)) / (interval_len + 1);  % tau[i - interval_len] + ... + tau[i]
    interval_len = fix(mean_tau / delt_t);  % tau不变的长度
    if interval_len + 5 > i
        continue
    end
    mean_tau = (sum_tau(i) - sum_tau(i - interval_len - 1)) / (interval_len + 1);  % tau[i - interval_len] + ... + tau[i]
    %var_tau = 1 / (r- l + 1)   *   sum_{i=l}^r  (a_i - a) ** 2
    %var_tau = 1 / (r- l + 1)   *   sum_{i=l}^r  (a_i**2 - 2a_i * a + a ** 2)
    %var_tau = 1 / (r- l + 1)   *  [ sum_{i=l}^r a_i ** 2      -  (r - l + 1) a ** 2 ]
    %var_tau = 1 / (r- l + 1)   *   sum_{i=l}^r a_i ** 2      -  a ** 2
    var_tau = (sum_tau2(i) - sum_tau2(i - interval_len - 1)) / (interval_len + 1) - mean_tau ^ 2;
    %var_tau = np.var(Dim4[-2, i - interval_len: i])
    pre_tau_pos = i - interval_len;
    if pre_tau_pos < 1
        continue
    end
    mean_len = fix(interval_len * 0.1);
    %now_mean = np.mean(Dim4[:-2, i - mean_len: i])
    %pre_mean = np.mean(Dim4[:-2, pre_tau_pos: pre_tau_pos + mean_len])
    now_mean = (sum_tau(i) - sum_tau(i - mean_len - 1)) / (mean_len + 1);  % tau[i - mean_len] + ... + tau[i]
    pre_mean = (sum_tau(pre_tau_pos + mean_len) - sum_tau(pre_tau_pos - 1)) / (mean_len + 1);  % tau[pre_tau_pos] + ... + tau[pre_tau_pos + mean_len]
    if abs(now_mean - pre_mean) > mean_base
        continue
    end
    if var_tau > var_base
        continue
    end
    now_state = Dim5(1:end-2, i);
    start_state = Dim5(1:end-2, pre_tau_pos);
    diff = sqrt(sum((now_state - start_state) .^ 2));
    if diff > diff_base
        continue
    end
    if mean_tau<1.5  %过滤掉小于1.5周期的假UPO
        continue
    end
    if isempty(window_TAU) || (min_diff > diff)
        min_diff = diff;
        window_TAU = [[mean_tau, mean_tau]];
        window_Time = [[t_now - mean_tau, t_now]];
        
        Diff = [diff];
        TAU=[mean_tau];
        Var=[var_tau];
        good = 1;
    end 
end
% TAU, Var, Diff, window_Time, window_TAU
end