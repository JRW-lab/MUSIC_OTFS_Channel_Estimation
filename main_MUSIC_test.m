clc; clear;

% Set paths and data
addpath(fullfile(pwd, 'Comm Functions'));
addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));

% Settings
init_frames = 10;
v_res = 0.01;
num_pilots_per_tsym = 1;
d_bs_to_road = 10; % in meters
                   % current method assumes travels towards BS for half the time
                   % then away from BS for the remaining time

% Set system parameters
new_frames = 100;
parameters = struct(...
    'system_name', "OTFS",...
    'CP', true,...
    'receiver_name', "SIC-MMSE",... 
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 20, ...
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

% Produce error if invalid number of pilot symbols
if mod(M/num_pilots_per_tsym,1) ~= 0 && M/num_pilots_per_tsym > 2*L+1
    error("Invalid number of pilot symbols, must divide N*M evenly.")
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
bit_errors = zeros(new_frames,1);
sym_errors = zeros(new_frames,1);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

% Generate delays and Dopplers
[~,chn_tau,chn_v] = channel_generation(Fc,vel);

% Parameters for change of Doppler w/ vehicle speed
t_frame = N*T;
t_travel = t_frame*new_frames;
d_traveled = 1000 * vel * t_travel / 3600; % in meters
starting_loc = -d_traveled / 2;
starting_angle = angle(starting_loc + 1j*d_bs_to_road);

% Parameters for channel estimation
block_size = M/num_pilots_per_tsym;
pilot_rows = zeros(num_pilots_per_tsym,1);
for i = 1:num_pilots_per_tsym
    pilot_rows(i) = block_size*i-L;
end
temp_mat = repmat(pilot_rows,[1,N]);
vec = (0:N-1);
temp_mat2 = temp_mat + vec*M;
pilot_idx = temp_mat2(:);
change_map = zeros(M,N);
for i = 1:num_pilots_per_tsym
    change_map((block_size*i-2*L):block_size*i,:) = 1;
end
change_map_vert = logical(change_map(:));

% Set variables for MUSIC
v_max = 100 * ceil((vel * (1000/3600)*Fc) / physconst('LightSpeed') / 100);
v_range = -v_max:v_res:v_max;
power_vec = pilot_idx-1;
e_temp = exp(1j*2*pi*Ts) .^ power_vec;
e = e_temp.^v_range;

% Sample frames for each diagonal entry
X_taps = zeros(N*num_pilots_per_tsym,init_frames,Lp-Ln+1);
H_cell = cell(init_frames,1);
nu_cell = cell(init_frames,1);
TX_bit_cell = cell(init_frames,1);
TX_sym_cell = cell(init_frames,1);
for frame = 1:new_frames

    if frame < init_frames
        fprintf("Generating frame %d/%d...\n",frame,new_frames)
    elseif frame == init_frames
        fprintf("Generating frame %d/%d and equalizing data until now...\n",frame,new_frames)
    else
        fprintf("Generating and equalizing frame %d/%d...\n",frame,new_frames)
    end

    % Generate data
    [TX_bit,TX_sym,x_DD] = gen_data(bit_order,S_alphabet,syms_per_f);

    % Add zeros to end of each pilot symbol block, insert pilot symbol per tsym
    TX_bit(change_map_vert,:) = -1;
    TX_sym(change_map_vert) = -1;
    x_DD(change_map_vert) = 0;
    x_DD(pilot_rows,:) = sqrt(N);
    if frame <= init_frames
        TX_bit_cell{frame} = TX_bit;
        TX_sym_cell{frame} = TX_sym;
    end

    % Shuffle data into needed formats
    % TX_bit_shuffled = Gamma_MN' * TX_bit;
    % TX_sym_shuffled = Gamma_MN' * TX_sym;
    % x_tilde = Gamma_MN' * x_DD;
    X_DD = reshape(x_DD,M,N);
    S = X_DD * F_N';
    s = S(:);

    % Shift Doppler according to vehicle motion
    ms_loc = starting_loc + t_frame * vel * (frame-1);
    new_angle = angle(ms_loc + 1j*d_bs_to_road);
    chn_v_sel = sin(new_angle) * chn_v / sin(starting_angle);

    % Generate channel
    t_offset = max_timing_offset * Ts;
    chn_g = channel_generation(Fc,vel);
    if shape == "rect" % rectangular ambiguity is closed form
        % Create H Matrix
        H = gen_H(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v_sel,shape,alpha,t_offset);
    else
        % Normalize tau and v to cohere with discrete ambig values
        chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
        chn_v_sel = round(chn_v_sel/res_chn_v)*res_chn_v;

        % Find direct tap indices and tap values
        l = (Ln:Lp).';
        tap_t_range = (l*Ts - chn_tau + t_offset) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_f_range = (ones(Lp-Ln+1,1) .* chn_v_sel) .* ones(Lp-Ln+1,length(chn_g),N*M);
        tap_t_range = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
        tap_f_range = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);
        tap_t_range(tap_t_range < 1) = 1;
        tap_t_range(tap_t_range > 1001) = 1001;
        tap_f_range(tap_f_range < 1) = 1;
        tap_f_range(tap_f_range > 1001) = 1001;

        % Create H Matrix
        H = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_v_sel,ambig_vals,tap_t_range,tap_f_range,t_offset);
    end

    % Makes all non-causal taps causal (easier for SIC-MMSE algorithm)
    H = circshift(H,-Ln);
    if frame <= init_frames
        H_cell{frame} = H;
    end

    % Generate shuffled DD noise, convert to time domain
    z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    z_DD = Gamma_MN * z_tilde;
    Z_DD = reshape(z_DD,M,N);
    Z_TF = F_M * Z_DD * F_N';
    W = F_M' * Z_TF;
    w = W(:);

    % Create receive vector
    nu = H * s;
    if frame <= init_frames
        nu_cell{frame} = nu;
    end

    % Collect samples from diagonals
    for i = 1:length(pilot_idx)
        x_samp = nu(pilot_idx(i):pilot_idx(i)+L);
        if frame <= init_frames
            X_taps(i,frame,:) = reshape(x_samp,1,1,length(x_samp));
        else
            X_taps(i,end,:) = reshape(x_samp,1,1,length(x_samp));
        end
    end

    % Begin estimation and equalization if initialization frames are complete
    if frame >= init_frames

        % Solve MUSIC for each diagonal
        P_scores = cell(Lp-Ln+1,1);
        v_est = cell(Lp-Ln+1,1);
        s_est = cell(Lp-Ln+1,1);
        X_cell = squeeze(num2cell(X_taps, [1 2]));
        if frame > init_frames
            X_taps = circshift(X_taps,[0,-1,0]);
        end
        for tap_idx = 1:Lp-Ln+1

            % Select current x samples
            X = X_cell{tap_idx};

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
            [scores,idx] = findpeaks(real(P_hat),'NPeaks',p_est,'SortStr','descend');
            P_scores{tap_idx} = scores.';
            v_est{tap_idx} = v_range(idx).';

            % "Re"create original Vandermonte matrix A (x = A * s), solve for s
            A = e_temp.^v_range(idx);
            s_est{tap_idx} = A \ X;

        end

        % Get combined set of Dopplers
        v_est_sorted = sort(unique(vertcat(v_est{:})));

        % Recreate A from Doppler estimates to get s vectors
        A_est = e_temp.^(v_est_sorted.');
        s_est = cellfun(@(x) lsqminnorm(A_est,x), X_cell, 'UniformOutput', false);

        % Create A again to estimate all taps
        power_vec_ext = (-Lp:(N*M-1-Ln))';
        e_temp2 = exp(1j*2*pi*Ts) .^ power_vec_ext;
        A_ext = e_temp2.^(v_est_sorted.');
        X_est = cellfun(@(x) A_ext * x, s_est, 'UniformOutput', false);

        % Channel reconstruction
        if frame == init_frames
            frame2_max = init_frames;
            frame2_min = 1;
        else
            frame2_max = frame;
            frame2_min = frame;
        end
        H_est = cell(init_frames,1);
        for frame2 = frame2_min:frame2_max

            % Select frame's RX vector
            if frame == init_frames
                nu_sel = nu_cell{frame2};
                TX_bit_sel = TX_bit_cell{frame2};
                TX_sym_sel = TX_sym_cell{frame2};
            else
                nu_sel = nu;
                TX_bit_sel = TX_bit;
                TX_sym_sel = TX_sym;
            end

            % Initialize H reconstruction and go by diagonal
            H_est_ext = zeros(N*M+L,N*M);
            for i = 1:length(X_est)

                % Select diagonal estimates
                X_sel = X_est{i};
                if frame == init_frames
                    X_diag = X_sel(:,frame2);
                else
                    X_diag = X_sel(:,end);
                end
                X_final = X_diag(Lp+1:(N*M+Lp));

                % Take beginning and end samples depending on causal or non-causal
                if i+Ln < 1 % Place non-causal part at beginning
                    X_end = X_diag(end+1+Ln:end+1-i);
                    X_final(1:(-Ln-i+1)) = X_end;
                else % Place causal part at end
                    X_end = X_diag(L-i-Ln:Lp);
                    X_final(end-i+Lp:end) = X_end;
                end

                % Reconstruct this diagonal, add to stack
                H_diag_mat = [zeros(i-1,N*M); diag(X_final); zeros(length(X_est)-i,N*M)];
                H_est_ext = H_est_ext + H_diag_mat;

            end

            % Make triangular part and make full H reconstruction
            H_end = H_est_ext((end-L+1):end,:);
            H_est = H_est_ext(1:N*M,:) + [H_end; zeros(N*M-L,N*M)];

            % Iterative Detector - JRW
            switch receiver_name
                case "SIC-MMSE"
                    [x_hat,iters_vec(frame2),t_RXiter_vec(frame2),t_RXfull_vec(frame2)] = equalizer_SIC_MMSE(nu_sel,H_est,N/num_pilots_per_tsym,M,2*L+1,Es,N0,S_alphabet,N_iters);
                otherwise
                    error("Unsupported receiver for the simulated system!")
            end

            % Hard detection for final x_hat
            dist = abs(x_hat.' - S_alphabet).^2;
            [~,min_index] = min(dist);
            RX_sym = min_index.' - 1;

            % Convert final RX_sym to RX_bit
            RX_bit = bit_order(RX_sym+1,:);

            % Error calculation
            if receiver_name == "SIC-MMSE"
                bit_error_vec = TX_bit_sel(~change_map_vert,:) ~= RX_bit(~change_map_vert,:);
                sym_error_vec = TX_sym_sel(~change_map_vert) ~= RX_sym(~change_map_vert);
            end
            bit_errors(frame2) = sum(bit_error_vec,"all");
            sym_errors(frame2) = sum(sym_error_vec,"all");
            if sum(bit_error_vec,"all") > 0
                frm_errors(frame2) = 1;
            else
                frm_errors(frame2) = 0;
            end

            % % DEBUGGING - Select current true H matrix
            % H_sel = H_cell{frame};
            % H_diag_true = diag(H_sel,-i+1);
            % norm(H_est{frame} - H_sel)

        end
    end
end

% Get parameters for throughput
frame_duration = N * T;
bandwidth_hz = M / T;

% Calculate BER, SER and FER
zero_syms = N*(2*L+1)*num_pilots_per_tsym;
metrics.BER = sum(bit_errors,"all") / (new_frames*(syms_per_f-zero_syms)*log2(M_ary));
metrics.SER = sum(sym_errors,"all") / (new_frames*(syms_per_f-zero_syms));
metrics.FER = sum(frm_errors,"all") / (new_frames);
metrics.Thr = (log2(M_ary) * (syms_per_f-zero_syms) * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
metrics.RX_iters = mean(iters_vec);
metrics.t_RXiter = mean(t_RXiter_vec);
metrics.t_RXfull = mean(t_RXfull_vec);