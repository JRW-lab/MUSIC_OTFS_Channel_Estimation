clc; clear;

% Set paths and data
addpath(fullfile(pwd, 'Comm Functions'));
addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));

% Settings
num_frames = 10;
sample_start = 1;
sample_rate = 6;
v_res = 0.01;
grouping_thresh = 8;
num_pilots_per_tsym = 1;

% Set system parameters
parameters = struct(...
    'system_name', "OTFS",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",... 
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 100000, ...
    'M', 64, ...
    'N', 16, ...
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
S_alphabet = sys.S;
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
Lp = sys.Lp;
Ln = sys.Ln;
if shape == "rect" || shape == "ideal"
    q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Make total channel length after delay shifting
L = Lp - Ln;

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
F_M = gen_DFT(M);

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
X_taps = zeros(N*num_pilots_per_tsym,num_frames,Lp-Ln+1);
for frame = 1:num_frames


    % Generate data
    [TX_bit,TX_sym,x_DD] = gen_data(bit_order,S_alphabet,syms_per_f);
    % % x_DD = [repmat([1; zeros(L,1)],floor(N*M/(L+1)),1); zeros(N*M-(L+1)*floor(N*M/(L+1)),1)];
    % % x_DD = zeros(N*M,1);
    % % x_DD(1:2:end) = 1;
    % if zero_padding
        % Count number of total zero padding needed
        zero_syms = N*L;

        % Add zeros to end of each time symbol
        change_map = zeros(M,N);
        change_map((M-2*L):M,:) = 1;
        % change_map_vert = logical(change_map(:));
        change_map_vert = true(M,N); % For now, just remove all data

        % Apply zeros over already-generated data
        TX_bit(change_map_vert,:) = -1;
        TX_sym(change_map_vert) = -1;
        x_DD(change_map_vert) = 0;

        % Insert pilot symbols L spacing apart per time symbol
        pilot_rows = zeros(num_pilots_per_tsym,1);
        for i = 1:num_pilots_per_tsym
            pilot_rows(i) = M-i*L-i+1;
        end
        % pilot_rows = [M-L,M-2*L-1,M-3*L-2];
        x_DD(pilot_rows,:) = sqrt(N);

    % else
    %     % Count number of total zero padding needed
    %     zero_syms = 0;
    % 
    %     % Create blank change map
    %     change_map_vert = false(M,N);
    % end
    TX_bit_shuffled = Gamma_MN' * TX_bit;
    TX_sym_shuffled = Gamma_MN' * TX_sym;
    x_tilde = Gamma_MN' * x_DD;
    X_DD = reshape(x_DD,M,N);
    S = X_DD * F_N';
    s = S(:);

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

    % Makes all non-causal taps causal (easier for SIC-MMSE algorithm)
    H = circshift(H,-Ln);

    % Generate shuffled DD noise, convert to time domain
    z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    z_DD = Gamma_MN * z_tilde;
    Z_DD = reshape(z_DD,M,N);
    Z_TF = F_M * Z_DD * F_N';
    W = F_M' * Z_TF;
    w = W(:);

    % Create receive vector
    nu = H * s;

    % Collect samples from diagonals
    pilot_idx = zeros(N,num_pilots_per_tsym);
    saved_taps = cell(N,num_pilots_per_tsym);
    for n = 0:N-1

        for i = 1:num_pilots_per_tsym
            pilot_idx(n+1,i) = (M*(n+1)-L*i-i+1);
            saved_taps{n+1,i} = nu(pilot_idx(n+1,i):pilot_idx(n+1,i)+L).';
        end
    end

    % Combine information for all pilot multiples per time sym
    pilot_idx_all = pilot_idx(:);
    saved_taps_all = saved_taps(:);
    saved_taps_mat = vertcat(saved_taps_all{:});

    % Sort all taps since I acquired them in reverse
    [~,sorted_idx] = sort(pilot_idx_all, 'ascend');
    pilot_idx_sorted = pilot_idx_all(sorted_idx);
    saved_taps_sorted = saved_taps_mat(sorted_idx,:);

    % Sample diagonals and take every Nth sample
    for l = 1:(Lp-Ln+1)
        x_samp = saved_taps_sorted(:,l);
        X_taps(:,frame,l) = x_samp;
    end

    % % Unwrap matrix to use all NM possible elements per diagonal
    % H_left = [zeros(Lp-Ln,N*M-(Lp-Ln)), H_shift(1:(Lp-Ln),N*M-(Lp-Ln)+1:N*M)];
    % H_center = H_shift;
    % H_center(1:(Lp-Ln),N*M-(Lp-Ln)+1:N*M) = 0;
    % H_unwrapped = [H_center; H_left];
    % 
    % % Sample diagonals and take every Nth sample
    % for l = 1:(Lp-Ln+1)
    %     h_diag = diag(H_unwrapped,-l+1);
    %     x_samp = h_diag(sample_start:sample_rate:N*M);
    %     X_taps(:,frame,l) = x_samp;
    % end

end

% Set variables for MUSIC
v_max = 100 * ceil((vel * (1000/3600)*Fc) / physconst('LightSpeed') / 100);
v_range = -v_max:v_res:v_max;
% power_vec = ((sample_start-1):sample_rate:(N*M-1))';
power_vec = pilot_idx_sorted-1;
e_temp = exp(1j*2*pi*Ts) .^ power_vec;
e = e_temp.^v_range;

% Solve MUSIC for each diagonal
P_scores = cell(Lp-Ln+1,1);
v_est = cell(Lp-Ln+1,1);
s_est = cell(Lp-Ln+1,1);
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

    % Estimate number of paths - Minimum Descriptive Length method
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
    Un = Us(:,(p_est+1):end);

    % Sweep through all possible Doppler values and generate costs
    norm_mat = Un' * e;
    d_sqrd = sum(abs(norm_mat).^2,1);
    P_hat = 1 ./ d_sqrd;

    % Estimate doppler shifts
    [P_scores{tap_idx},idx] = findpeaks(real(P_hat),'NPeaks',p_est,'SortStr','descend');
    v_est{tap_idx} = v_range(idx).';

    % "Re"create original Vandermonte matrix A (x = A * s), solve for s
    A = e_temp.^v_range(idx);
    s_est{tap_idx} = A \ X;

end

% A_tap_range = (0:(N*M-1)).';
A_tap_range = (sample_start:sample_rate:(N*M-1)).';

% 
thresh = 1;
v_est_true = [];
for i = 1:length(v_est)
    v_est_sel = v_est{i};

    if isempty(v_est_true) || isempty(v_est_sel)
        v_est_true = [v_est_true; v_est_sel];
    else
        dist = (v_est_true - v_est_sel.');
        min_dist = min(abs(dist).^2,[],1);
        new_locs = min_dist > thresh;

        v_est_true = [v_est_true; v_est_sel(new_locs)];
    end
end
v_est_sorted = sort(v_est_true);


%%

% Collect all estimated Dopplers
v_est_all = vertcat(v_est{:});
v_est_all = sort(v_est_all);

% Group by similarity and only take unique estimates
groups = [true; abs(diff(v_est_all)) > grouping_thresh];
v_est_sorted = v_est_all(groups);
chn_v_sorted = sort(chn_v).';


% Just print the values so we can see them for now
v_est_sorted.'
chn_v_sorted.'