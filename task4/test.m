function run_ocfe_collocation_demo()
%% ============================================================
%  OCFE (Radau, degree=3) 全离散化 + 顺序单元求解 + 与 ode15s 对比
%  过程：A <- B -> C（按你给的动力学），目标：高 cA
%  时域：0 ~ 3600 s；控制每 600 s 一段（共 6 段），默认恒 343 K
% ============================================================

clf; close all;

%% ------------------- 基本设置 -------------------
TF      = 3600;            % 1小时（秒）
DEG     = 3;               % collocation 阶数（Radau根，三阶）
NE_CAND = [6, 12, 24, 36, 60, 90];   % 候选有限元个数（你可改）
TSTEP   = 600;             % 控制分段长度（秒），10 分钟
NBLK    = TF / TSTEP;      % 控制段数，题目=6

% 初始条件（题给）
y0 = [10; 18; 0; 343];     % [cA, cB, cC, T]
T0_init = 343; Tj_init = 343;

% 分段常数控制（示例：整段恒定 343 K；你可改成 6 段向量）
uT0 = T0_init*ones(NBLK,1);
uTj = Tj_init*ones(NBLK,1);

% 控制边界 & 步进限制（供你后续做 OCP 用；本示例先只用于检查）
U_MIN = 273; U_MAX = 373;
MAX_STEP = 5; % K

%% ------------------- 基准解（常规 ODE） -------------------
% 为了比较 OCFE 的精度，我们用 ode15s 做一个“参考解”
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
sol  = ode15s(@(t,y) reactor_ode(t,y,uT0,uTj,TSTEP), [0 TF], y0, opts);

% 用等间隔网格采样做对比（每分钟）
t_grid = linspace(0,TF,61);
y_ref  = deval(sol, t_grid);

%% ------------------- 误差阈值与自动选取 NE -------------------
% 误差阈值（可按需要调）：cA,cB,cC 用 1e-2；温度用 0.2 K
tol_x  = [1e-2; 1e-2; 1e-2; 2e-1];

chosen_NE = NE_CAND(end); % 默认用最后一个（若没人达标）
for NE = NE_CAND
    [t_col, y_col, Tmax] = simulate_by_ocfe(NE, DEG, TF, y0, uT0, uTj, TSTEP);
    % 插值到 t_grid 对比
    y_col_grid = interp1(t_col, y_col.', t_grid, 'pchip').';

    err = max(abs(y_col_grid - y_ref), [], 2); % 各状态最大绝对误差
    if all(err <= tol_x)
        chosen_NE = NE;
        fprintf('选择 NE = %d（满足误差阈值）\n', NE);
        break;
    else
        fprintf('NE = %d 未通过：误差 = [%g, %g, %g, %g]\n', NE, err);
    end
    if Tmax > 350
        fprintf(2,'  警告：NE=%d 下出现 T=%.2f K 超过 350K（仅提示）\n', NE, Tmax);
    end
end

%% ------------------- 用选定的 NE 重算并绘图 -------------------
[t_col, y_col, Tmax] = simulate_by_ocfe(chosen_NE, DEG, TF, y0, uT0, uTj, TSTEP);
y_col_grid = interp1(t_col, y_col.', t_grid, 'pchip').';

figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on;
plot(t_grid/60, y_ref(1,:), '-', 'LineWidth',1.5);
plot(t_grid/60, y_col_grid(1,:), '--', 'LineWidth',1.5);
ylabel('c_A (mol/m^3)'); xlabel('Time (min)'); legend('ode15s','OCFE'); title('c_A');

nexttile; hold on; box on;
plot(t_grid/60, y_ref(2,:), '-', 'LineWidth',1.5);
plot(t_grid/60, y_col_grid(2,:), '--', 'LineWidth',1.5);
ylabel('c_B (mol/m^3)'); xlabel('Time (min)'); legend('ode15s','OCFE'); title('c_B');

nexttile; hold on; box on;
plot(t_grid/60, y_ref(3,:), '-', 'LineWidth',1.5);
plot(t_grid/60, y_col_grid(3,:), '--', 'LineWidth',1.5);
ylabel('c_C (mol/m^3)'); xlabel('Time (min)'); legend('ode15s','OCFE'); title('c_C');

nexttile; hold on; box on;
plot(t_grid/60, y_ref(4,:), '-', 'LineWidth',1.5);
plot(t_grid/60, y_col_grid(4,:), '--', 'LineWidth',1.5);
yline(350,'k:','350K limit');
ylabel('T (K)'); xlabel('Time (min)'); legend('ode15s','OCFE'); title(sprintf('T (NE=%d)', chosen_NE));

sgtitle('OCFE（Radau-3）与 ODE 仿真对比');

if Tmax > 350
    fprintf(2,'最终选定 NE=%d 下最大过程温度 %.2f K 超过 350 K（仅提示）。\n', chosen_NE, Tmax);
end

%% -------------------（可选）检查控制步进与边界 -------------------
chk_ctrl('T0', uT0, U_MIN, U_MAX, MAX_STEP);
chk_ctrl('Tj', uTj, U_MIN, U_MAX, MAX_STEP);

end % ====== 主函数结束 ======



%% ========= 工具函数：逐单元 OCFE 仿真（Radau-III, degree=3） =========
function [t_all, y_all, Tmax] = simulate_by_ocfe(NE, DEG, TF, y0, uT0, uTj, TSTEP)
    % Collocation 系数（Radau 根，含 0 端点）
    [Cmat, Dvec, tau_root] = radau_coeffs(DEG); % C: d x (d+1), D: (d+1) x 1
    h  = TF / NE;                               % 每单元时长
    nx = 4; d = DEG;

    % 结果存储：把每个单元的 collocation 点与末端都存下来
    t_all = 0;
    y_all = y0(:).';
    yk    = y0(:);   % 当前单元起点状态
    Tmax  = y0(4);

    opts_fsolve = optimoptions('fsolve', ...
        'Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxIterations',400);

    for k = 1:NE
        t0 = (k-1)*h;

        % 本单元的 collocation 点时间（相对全局）
        t_c = t0 + h * tau_root(2:end); % 不包含 0，只要 3 个内部点（含1）
        % 控制段索引（按点映射）；末端点也用同段
        blk_idx = min(numel(uT0), floor(t_c / TSTEP) + 1);
        blk_end = min(numel(uT0), floor((t0+h)/TSTEP) + 1);

        % 未知向量：Z = [z1; z2; z3; x_end]，每段 z_j 与 x_end 都是 4 维
        Z0 = repmat(yk, d+1, 1); % 初值猜测：全设为起点
        fun = @(Z) residual_element(Z, yk, h, Cmat, Dvec, nx, d, ...
                                    uT0(blk_idx), uTj(blk_idx), uT0(blk_end), uTj(blk_end));

        Zsol = fsolve(fun, Z0, opts_fsolve);
        % 取出 z_j 与 x_{k+1}
        Zstg  = reshape(Zsol(1:d*nx), nx, d);
        ykp1  = Zsol(d*nx+1:end);

        % 更新轨迹
        t_all = [t_all; t_c(:); t0+h];
        y_all = [y_all; Zstg.'; ykp1(:).'];
        yk    = ykp1(:);
        Tmax  = max(Tmax, max(Zstg(4,:), [], 'all'));
        Tmax  = max(Tmax, yk(4));
    end
end

function R = residual_element(Z, yk, h, C, D, nx, d, uT0_blk, uTj_blk, uT0_end, uTj_end)
    % 把未知量拆成各 stage 与末端
    Z = Z(:);
    Zstg = reshape(Z(1:d*nx), nx, d); % [x1..x4]x3
    xend = Z(d*nx+1:end);

    % 组装“节点+stage”便于使用 Lagrange 系数组合：X0,X1,X2,X3
    X = [yk, Zstg];  % size nx x (d+1)

    % 1) collocation 方程：对每个 j=1..d： sum_r C(j,r)*X_r - h*f(Z_j, u) = 0
    Rcol = zeros(nx*d,1);
    for j = 1:d
        % 组合项 sum_r C(j,r)*X_r
        comb = zeros(nx,1);
        for r = 1:(d+1)
            comb = comb + C(j,r) * X(:,r);
        end
        f = reactor_rhs(Zstg(:,j), uT0_blk(j), uTj_blk(j)); % 该 collocation 点用对应控制段
        Rcol((j-1)*nx+(1:nx)) = comb - h * f;
    end

    % 2) 末端连续性： sum_r D(r)*X_r - x_{k+1} = 0
    xcont = zeros(nx,1);
    for r = 1:(d+1)
        xcont = xcont + D(r) * X(:,r);
    end
    Rend = xcont - xend;

    R = [Rcol; Rend];
end

%% ========= Radau-III（degree=3）Lagrange 基下的 C、D 系数 =========
function [C, D, tau_root] = radau_coeffs(d)
    % 只实现 d=3；若你要通用化，可按 Lagrange 基的通式生成
    if d ~= 3
        error('此示例实现的是三阶 Radau（d=3）。');
    end
    % Radau IIA 节点（含 1.0），再加 0 作为 tau_0（与文献/经典实现一致）
    tau = [0.155051025721680, 0.644948974278317, 1.0]; % 3 个内部点（含末端 1）
    tau_root = [0, tau];                                % 含起点 0

    % 构造 Lagrange 多项式，计算 C(j,r)=L_r'(tau_j)、D(r)=L_r(1)
    d1 = d + 1;
    C  = zeros(d, d1);
    D  = zeros(d1,1);

    % 逐个基函数 L_r
    for r = 1:d1
        % L_r(τ) = Π_{m≠r} (τ - τ_m)/(τ_r - τ_m)
        % L_r'(τ) = Σ_{m≠r} [ 1/(τ_r-τ_m) * Π_{n≠r,n≠m} (τ - τ_n)/(τ_r - τ_n) ]
        % 数值评估函数：
        L  = @(t) lag_basis_val(t, r, tau_root);
        Ld = @(t) lag_basis_der(t, r, tau_root);

        % D(r) = L_r(1)
        D(r) = L(1.0);

        % C(j,r) = L_r'(tau_j)  （j=1..d；对应 tau_root(2:end)）
        for j = 1:d
            C(j,r) = Ld(tau_root(j+1));
        end
    end
end

function val = lag_basis_val(t, r, tau_root)
    d1  = numel(tau_root);
    tr  = tau_root(r);
    val = 1.0;
    for m = 1:d1
        if m == r, continue; end
        val = val * (t - tau_root(m)) / (tr - tau_root(m));
    end
end

function val = lag_basis_der(t, r, tau_root)
    % L_r'(t) 的直接公式（按乘积求导展开）
    d1  = numel(tau_root);
    tr  = tau_root(r);
    val = 0.0;
    for m = 1:d1
        if m == r, continue; end
        term = 1.0 / (tr - tau_root(m));
        for n = 1:d1
            if n == r || n == m, continue; end
            term = term * (t - tau_root(n)) / (tr - tau_root(n));
        end
        val = val + term;
    end
end

%% ========= 反应器 ODE（按题）供 ode15s 使用：u 分段常数 =========
function dydt = reactor_ode(t, y, uT0, uTj, TSTEP)
    % 找到当前控制段
    idx = min(numel(uT0), floor(t / TSTEP) + 1);
    T0  = uT0(idx);
    Tj  = uTj(idx);
    dydt = reactor_rhs(y, T0, Tj);
end

%% ========= 反应器 RHS：与题中 daeSystemLHS 等价（去掉 X 和 PARAM 的壳）=========
function dydt = reactor_rhs(y, T0, Tj)
    % 常数
    EA      = 69000.0;    % J/mol
    EB      = 72000.0;    % J/mol
    kA0     = 5.0e6;      % 1/s
    kB0     = 1.0e7;      % 1/s
    rho     = 800.0;      % kg/m^3
    cp      = 3.5;        % J/kg/K 
    UA      = 1.4;        % W/K
    deltaHA =  45000.0;   % J/mol   
    deltaHB = -55000.0;   % J/mol
    F       = 6.5e-4;     % m^3/s
    V       = 1.0;        % m^3
    CA0     = 5.0;        % mol/m^3
    CB0     = 15.0;       % mol/m^3

    cA = y(1); cB = y(2); cC = y(3); T = y(4);

    kA = kA0 * exp(-EA / (8.314 * T));
    kB = kB0 * exp(-EB / (8.314 * T));

    dcAdt = (F/V)*(CA0 - cA) - kA*cA + kB*cB;
    dcBdt = (F/V)*(CB0 - cB) - kB*cB; 
    dcCdt = -(F/V)*cC + kA*cA;

    dTdt  = (F/V)*(T0 - T) ...
          + (UA / (rho * cp * V)) * (Tj - T) ...
          + (-deltaHA * kA * cA + (-deltaHB) * kB * cB) / (rho * cp);

    dydt = [dcAdt; dcBdt; dcCdt; dTdt];
end

%% ========= 小工具：检查控制边界与步进限制 =========
function chk_ctrl(name, u, umin, umax, maxstep)
    if any(u < umin) || any(u > umax)
        fprintf(2,'%s 有越界值（应在 [%g,%g]）\n', name, umin, umax);
    end
    if any(abs(diff(u)) > maxstep)
        fprintf(2,'%s 某步进超过 %g K 限制\n', name, maxstep);
    end
end
