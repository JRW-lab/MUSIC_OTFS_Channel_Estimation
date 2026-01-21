clc; clear;

% Set paths and data
addpath(fullfile(pwd, 'Comm Functions'));
addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));

% Settings
num_frames = 10;
sample_start = 4;
sample_rate = 9;
v_res = 0.01;
grouping_thresh = 5;

% Set system parameters
parameters = struct(...
    'system_name', "OTFS",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",... 
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 100, ...
    'M', 32, ...
    'N', 32, ...
    'T', 1 / 15000, ...
    'Fc', 4e9, ...
    'vel', 120, ...
    'shape', "rrc", ...
    'alpha', 0.4, ...
    'Q', 2);

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Change settings if CP
if ~CP
    error("This method was derived with CP implicitly encoded in the channel structure.")
end

% Input conversion
sys = obj_comms_OTFS;
sys.EbN0_db = EbN0;
sys.M_ary = M_ary;
sys.N_tsyms = N;
sys.M_sbcars= M;
sys.filter = shape;
sys.q = Q;
sys.sbcar_spacing = 1 / T;
sys.Fc = Fc;
sys.v_vel = vel;
sys.rolloff = alpha;

% Equalizer settings
N_iters = 8;
ambig_res = 101;

% Inputs from system object
Es = sys.Es;
S = sys.S;
M_ary = sys.M_ary;
N0 = sys.N0;
M = sys.M_sbcars;
N = sys.N_tsyms;
Fc = sys.Fc;
vel = sys.v_vel;
if M_ary == 2
    bit_order = [0;1];
elseif M_ary == 4
    bit_order = [0,0;0,1;1,0;1,1];
end
shape = sys.filter;
alpha = sys.rolloff;
q = sys.q;
T = sys.T;
Ts = sys.Ts;
Lp = sys.Lp + 1;
Ln = sys.Ln - 1;
if shape == "rect" || shape == "ideal"
    q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

load_test = sys.ambig_vals;
if isempty(load_test)
    if shape ~= "rect"
        if ~exist("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables", 'dir')
            mkdir("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables")
        end
        filename = sprintf("Pre-rendered Lookup Tables\\OTFS-DD Cross-Ambiguity Tables\\ambig_discrete_M%d_N%d_T%d_Fc%d_vel%d_%s_alpha%.1f_q%d.mat",M,N,T,Fc,vel,shape,alpha,q);
        if isfile(filename)
            loaded_file = load(filename);
            ambig_t_range = loaded_file.ambig_t_range;
            ambig_f_range = loaded_file.ambig_f_range;
            ambig_vals = loaded_file.ambig_vals;
        else
            ambig_t_lim = Ts*(q + ceil((2510*10^(-9))/Ts));
            ambig_f_lim = ((vel * (1000/3600))*Fc) / (physconst('LightSpeed'));
            ambig_t_range = linspace(-ambig_t_lim,ambig_t_lim,ambig_res);
            ambig_f_range = linspace(-ambig_f_lim,ambig_f_lim,ambig_res);
            ambig_vals = zeros(ambig_res);
            % update_vals = floor(length(ambig_t_range)*(linspace(.01,.99,34)));
            for k = 1:length(ambig_t_range)
                for l = 1:length(ambig_f_range)
                    ambig_vals(k,l) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,shape,alpha,q,ambig_res);
                end
            end
            % fprintf("\n")

            % Save to file
            save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
        end

        % fprintf("Complete!\n\n")
    else
        ambig_t_range = [];
        ambig_f_range = [];
        ambig_vals = [];
    end
    sys.ambig_t_range = ambig_t_range;
    sys.ambig_f_range = ambig_f_range;
    sys.ambig_vals = ambig_vals;
end

load_test = sys.R;
if isempty(load_test)
    % Set up noise covariance matrix if not done already
    if shape ~= "rect"
        if ~exist("Pre-rendered Lookup Tables\\OTFS-DD Noise Covariance Matrices", 'dir')
            mkdir("Pre-rendered Lookup Tables\\OTFS-DD Noise Covariance Matrices")
        end
        filename = sprintf("Pre-rendered Lookup Tables\\OTFS-DD Noise Covariance Matrices\\Rzddt_T%d_N%d_M%d_q%d_%s_alpha%.1f.mat",T,N,M,q,shape,alpha);
        if isfile(filename)
            loaded_file = load(filename);
            R_half = loaded_file.R_half;
            R = loaded_file.R;
        else
            % fprintf("Generating noise covariance matrix (one time process)...\n")
            [R_half,R] = gen_Rzddt(T,N,M,q,ambig_res,shape,alpha);
            % fprintf("\n")

            % Save to file
            save(filename,"R_half","R");
        end
    else
        R_half = eye(N*M);
        R = R_half;
    end
    sys.R_half = R_half;
    sys.R = R;
end

% Needed variables and matrices setup
syms_per_f = M*N;
Gamma_MN = gen_Gamma_MN(M,N);
F_N = gen_DFT(N);

% Bring in previously loaded data
ambig_t_range = sys.ambig_t_range;
ambig_f_range = sys.ambig_f_range;
ambig_vals = sys.ambig_vals;
R = sys.R;
R_half = sys.R_half;
if shape ~= "rect"
    res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
    res_chn_v = ambig_f_range(2)-ambig_f_range(1);
end

% Reset bit errors for each SNR
bit_errors = zeros(num_frames,syms_per_f*log2(M_ary));
sym_errors = zeros(num_frames,syms_per_f);
frm_errors = zeros(num_frames,1);
iters_vec = zeros(num_frames,1);
t_RXiter_vec = zeros(num_frames,1);
t_RXfull_vec = zeros(num_frames,1);

% Generate delays and Dopplers
[~,chn_tau,chn_v] = channel_generation(Fc,vel);

% Sample frames for each diagonal entry
X_taps = zeros(length(sample_start:sample_rate:N*M),num_frames,Lp-Ln+1);
for frame = 1:num_frames

    % Generate channel
    t_offset = max_timing_offset * Ts;
    chn_g = channel_generation(Fc,vel);
    if shape == "rect" % rectangular ambiguity is closed form
        % Create H Matrix
        H = gen_H(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha,t_offset);
    else
        % Normalize tau and v to cohere with discrete ambig values
        chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
        chn_v = round(chn_v/res_chn_v)*res_chn_v;

        % Find direct tap indices and tap values
        l = (Ln:Lp).';
        tap_t_range = (l*Ts - chn_tau + t_offset) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_t_range = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
        tap_f_range = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);
        tap_t_range(tap_t_range < 1) = 1;
        tap_t_range(tap_t_range > 1001) = 1001;
        tap_f_range(tap_f_range < 1) = 1;
        tap_f_range(tap_f_range > 1001) = 1001;

        % Create H Matrix
        H = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_v,ambig_vals,tap_t_range,tap_f_range,t_offset);
    end

    % Unwrap matrix to use all NM possible elements per diagonal
    H_shift = circshift(H,-Ln); % Makes all non-causal taps causal
    H_left = [zeros(Lp-Ln,N*M-(Lp-Ln)), H_shift(1:(Lp-Ln),N*M-(Lp-Ln)+1:N*M)];
    H_center = H_shift;
    H_center(1:(Lp-Ln),N*M-(Lp-Ln)+1:N*M) = 0;
    H_unwrapped = [H_center; H_left];

    % Sample diagonals and take every Nth sample
    for l = 1:(Lp-Ln+1)
        h_diag = diag(H_unwrapped,-l+1);
        x_samp = h_diag(sample_start:sample_rate:N*M);
        X_taps(:,frame,l) = x_samp;
    end

end

% Set variables for MUSIC
v_max = 100 * ceil((vel * (1000/3600)*Fc) / physconst('LightSpeed') / 100);

% Solve MUSIC for each diagonal
v_est = cell(Lp-Ln+1,1);
for tap_idx = 1:Lp-Ln+1

    % Select current x samples
    X = X_taps(:,:,tap_idx);

    % Skip this diagonal if it doesn't have sufficient information
    if norm(X) < 1e-15
        continue;
    end

    % Generate covariance and select noise subspace
    R_x = X * X' / size(X,2);
    [U,D] = eig(R_x);
    [~, ind] = sort(diag(D), 'descend');
    Us = U(:,ind);

    % Estimate number of paths - Minimum Descriptive Length method (UNFINISHED)
    p = size(D,1);
    N_x = size(X,2);
    MDL = zeros(p,1);
    eigvals = abs(diag(D));
    for k = 1:p
        % Find products for numerator and denominator
        prod_num = 0;
        prod_den = 0;
        for i = k+1:p
            prod_num = prod_num + eigvals(i)^(1/(p-k));
            prod_den = prod_den + eigvals(i);
        end
        prod_den = prod_den / (p-k);

        % Find MDL(k) metric
        MDL(k) = -(p-k)*N_x*log(prod_num/prod_den) + 0.5*k*(2*p-k)*log(N);
    end
    [~,p_est] = min(abs(MDL));
    % p_est = length(chn_v);
    Un = Us(:,(p_est+1):end);

    % Sweep through all possible Doppler values and generate costs
    v_range = -v_max:v_res:v_max;
    power_vec = ((sample_start-1):sample_rate:(N*M-1))';
    e_temp = exp(1j*2*pi*Ts) .^ power_vec;
    e = e_temp.^v_range;
    norm_mat = Un' * e;
    d_sqrd = sum(abs(norm_mat).^2,1);
    P_hat = 1 ./ d_sqrd;

    % % Plot cost function (-P_hat)
    % figure(1)
    % semilogy(v_range,-P_hat)
    % xlabel("Estimated Doppler")
    % ylabel("Cost function")
    % 
    % 1;

    % Estimate doppler shifts
    [~,idx] = findpeaks(real(P_hat),'NPeaks',p_est,'SortStr','descend');
    v_est{tap_idx} = v_range(idx).';

    1;

end

% Collect all estimated Dopplers
v_est_all = vertcat(v_est{:});
v_est_all = sort(v_est_all);

% Group by similarity and only take unique estimates
groups = [true; abs(diff(v_est_all)) > grouping_thresh];
v_est_sorted = v_est_all(groups);

% Just print the values so we can see them for now
v_est_sorted.'
sort(chn_v)

1;
