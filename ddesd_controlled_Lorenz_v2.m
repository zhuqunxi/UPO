% 具有状态依赖时滞的 DDE ddesd  https://ww2.mathworks.cn/help/matlab/math/state-dependent-delay-problem.html
% 参考文献：Leveraging neural differential equations and adaptive delayed feedback to detect unstable periodic orbits based on irregularly sampled time series，https://aip.scitation.org/doi/full/10.1063/5.0143839
%  suggestion: INIT_GAMA=0.01, R1=R2=0.1
%Example: ddesd_controlled_Lorenz_v2(0, 0, 1, 1, 0.01, 0.1, 0.1)
function ddesd_controlled_Lorenz_v2(node, init, control, init_tau, init_gama, r1, r2) 
global NODE INIT CONTROL INIT_TAU INIT_GAMA R1 R2
NODE = node; % 0: Trus system, 1: Neural ODE
INIT = init; % 0: 不用历史数据， 1：使用历史数据
CONTROL = control; % 0: 不受控，跑真实系统， 1：受控，跑受控系统
INIT_TAU = init_tau; % ADFC 可调参数
INIT_GAMA = init_gama;  % ADFC 可调参数
R1 = r1;  % ADFC 可调参数
R2 = r2;  % ADFC 可调参数
diff_base = 0.5; % 用来数值上判定是否为UPO参数，包括下面的参数
var_base = 1e-5;
mean_base = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Initial_T Initial_X
N_init = 1000;
Init_data_file = './Prediction/NODE & UPO/UPO_long_1216/dataset/data_lorenz/';
Initial_T = readtable(sprintf('%s%s', Init_data_file, 'time_point.csv'));
Initial_X = readtable(sprintf('%s%s', Init_data_file, 'trajectory.csv'));
% Initial_T = Initial_T{2:N_init, 2} - Initial_T{N_init, 2};  % t in [-a, 0], a>TAU
% Initial_X = Initial_X{2:N_init, 2};
Initial_T = Initial_T{end - N_init:end, 2} - Initial_T{end, 2};  % t in [-a, 0], a>TAU
Initial_X = Initial_X{end - N_init:end, 2:4};
interp1(Initial_T, Initial_X, -0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%RelTol = 1e-4; % 1e-3
%AbsTol = 1e-8; % 1e-6
RelTol = 1e-3; % Solver的误差控制参数
AbsTol = 1e-6;
options = ddeset('RelTol', RelTol, 'AbsTol', AbsTol);
Time_Scale = 100;
%tspan = [0 50];
T = 310;
tspan = linspace(0, T,  T * Time_Scale + 1);

if node
    global model_state_dict
    model_state_dict = load('./model_mat/model_state_dict_Lorenz.mat'); %训练好的神经网络
    model_state_dict.input_layer_0_weight = double(model_state_dict.input_layer_0_weight);
    model_state_dict.input_layer_0_bias = double(model_state_dict.input_layer_0_bias);
    model_state_dict.input_layer_2_weight = double(model_state_dict.input_layer_2_weight);
    model_state_dict.input_layer_2_bias = double(model_state_dict.input_layer_2_bias);
    model_state_dict.output_layer_0_weight = double(model_state_dict.output_layer_0_weight);
    model_state_dict.output_layer_0_bias = double(model_state_dict.output_layer_0_bias);
    model_state_dict.output_layer_2_weight = double(model_state_dict.output_layer_2_weight);
    model_state_dict.output_layer_2_bias = double(model_state_dict.output_layer_2_bias);
    model_state_dict;
end
sol = ddesd(@ddefun, @dely, @history, tspan, options);

% 
yint = deval(sol,tspan);

mask =  sprintf('control_%d_node_%d_init_%d_itau_%s_igama%s_r1_%s_r2_%s',CONTROL, NODE, INIT, num2str(INIT_TAU), num2str(INIT_GAMA), num2str(R1), num2str(R2));
    
global TAU Var Diff window_Time window_TAU
good = Find_Period(tspan, yint, diff_base, var_base, mean_base);
fprintf('\ngood %d  %s\n', good, mask)
if good %找到UPO
    fprintf('********************************\n')
    fprintf('Tau: %f\n', TAU(1)) % UPO的周期
    fprintf('Var: %e\n', Var(1)) % 判定UPO稳定的方差
    fprintf('State Difference: %f\n', Diff(1)) %判定UPO首尾的状态误差
    fprintf('********************************\n')
    delt_t = tspan(2) - tspan(1); 
    index_st = fix(window_Time(1) / delt_t) + 1;
    index_ed = fix(window_Time(2) / delt_t) + 1;
    
    
    fig = figure;
    subplot(2,1,1)
    if CONTROL
        plot3(yint(1,index_st:index_ed),yint(2,index_st:index_ed),yint(3,index_st:index_ed), '-')
    else
        plot3(yint(1, :),yint(2, :),yint(3, :), '-')
    end

    
    % plot3(yint(1,:),yint(2,:),yint(3,:), '-')
    title(mask)
    % axis([0, 1.5, 0, 1.5])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view([45,15])
    grid on

    subplot(2,1,2)
    hold on
    plot(tspan, yint(end - 1, :))
    xlabel('t')
    ylabel('\tau(t)')
    plot(window_Time, window_TAU, 'r-')
    title(sprintf('Tau: %f, Var: %e, State Diff: %f', TAU(1), Var(1), Diff(1)))
    
    saveas(fig, sprintf('./data/Lorenz/node_%d/fig/%s.jpg', node, mask))
    close(fig)
    data_file = sprintf('./data/Lorenz/node_%d/%s.mat', node, mask);
    save(data_file,'tspan', 'yint', 'CONTROL', 'NODE', 'INIT', 'INIT_TAU', 'INIT_GAMA', 'R1', 'R2', 'TAU', 'Var', 'Diff', 'window_Time', 'window_TAU')

end


end

% 编写时滞代码
% 首先，编写一个函数来定义方程组中的时滞。此方程组中唯一存在的时滞位于项 ?y2{e^[1?y2(t)]} 中
function d = dely(t,y)
d = [t - y(end-1)];
end


% 求解器会自动将这些输入传递给该函数，但是变量名称决定如何编写方程代码。在这种情况下
% Z(2,1) → y2{e^[1?y2(t)]}
function dydt = ddefun(t,y,Z)
global NODE CONTROL R1 R2

r1 = R1;
r2 = R2;


diffs = Z(2, 1) - y(2);
if NODE
    dxdt = odefun_node(y(1:3));
else
    dxdt = odefun_true(y(1:3));
end

C_0 = 0.1;
S_t = y(end) * diffs;
if S_t > C_0
    S_t = C_0;
end
if S_t < -C_0
    S_t = -C_0;
end
    
if CONTROL
%     S_t = 0; % 保证非侵入性  控制很小
    dxdt(2) = dxdt(2) + S_t;
%     dxdt(2) = dxdt(2) + y(end) * diffs;
end


dtaudt = -r1 * diffs;
dgamadt = r2 * diffs ^ 2;
dydt = [dxdt(1);dxdt(2);dxdt(3);dtaudt; dgamadt];
end

function dydt = odefun_node(y)
global model_state_dict
dx = y'; % [x(t), x(t-tau)];   size(y) = [3, 1]  需要转置
% size(dx)
% size(model_state_dict.input_layer_0_weight)
% size(model_state_dict.input_layer_0_bias)

dx = dx * model_state_dict.input_layer_0_weight' + model_state_dict.input_layer_0_bias;
dx = tanh(dx);
dx = dx * model_state_dict.input_layer_2_weight' + model_state_dict.input_layer_2_bias;
dx = tanh(dx);
dx = dx * model_state_dict.output_layer_0_weight' + model_state_dict.output_layer_0_bias;
dx = tanh(dx);
dx = dx * model_state_dict.output_layer_2_weight' + model_state_dict.output_layer_2_bias;

dydt = dx';
end

function dydt = odefun_true(y)
rho = 28.0;
sigma = 10.0;
beta = 8.0 / 3.0;
dydt = [sigma * (y(2) - y(1));  - sigma * y(1) * y(3) + rho * y(1) - y(2); sigma * y(1) * y(2) - beta * y(3)];
end


% 编写历史解代码
% 接下来，创建一个函数来定义历史解。历史解是时间 t≤t0 的解。
function v = history(t) % history function for t < t0
% global Init_tau
% v = [1.254204;	1.46407;	1.0706747; 1; 0.01];
% v = [1.4003375;	1.2110717;	3.6193776; Init_tau; 0.01];

global INIT Initial_T Initial_X INIT_TAU INIT_GAMA
if INIT
    X = interp1(Initial_T, Initial_X, t);
    v = [X(1); X(2); X(3); INIT_TAU; INIT_GAMA]; % INIT_GAMA 0.01 suggested
else
    v = [1.4003375;	1.2110717;	3.6193776; INIT_TAU; INIT_GAMA];
end

end

