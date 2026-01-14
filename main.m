clc; clear;

% Set paths and data
% addpath(fullfile(pwd, 'Meta Functions'));
addpath(fullfile(pwd, 'Comm Functions'));
addpath(fullfile(pwd, 'Comm Functions/Custom Functions'));
addpath(fullfile(pwd, 'Comm Functions/Generation Functions'));
% addpath(fullfile(pwd, 'Comm Functions/OFDM Functions'));
% addpath(fullfile(pwd, 'Comm Functions/OTFS Functions'));
addpath(fullfile(pwd, 'Comm Functions/OTFS-DD Functions'));
% addpath(fullfile(pwd, 'Comm Functions/ODDM Functions'));
% addpath(fullfile(pwd, 'Comm Functions/TODDM Functions'));
addpath(fullfile(pwd, 'Comm Functions/TX RX Functions'));

% Set parameters and number of frames
new_frames = 1000;
parameters = struct(...
    'system_name', "OTFS",...
    'CP', true,...
    'receiver_name', "CMC-MMSE",... 
    'max_timing_offset', 0.0,...
    'M_ary', 4, ...
    'EbN0', 100, ...
    'M', 8, ...
    'N', 3, ...
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
% rng('default');

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
bit_errors = zeros(new_frames,syms_per_f*log2(M_ary));
sym_errors = zeros(new_frames,syms_per_f);
frm_errors = zeros(new_frames,1);
iters_vec = zeros(new_frames,1);
t_RXiter_vec = zeros(new_frames,1);
t_RXfull_vec = zeros(new_frames,1);

% Simulation loop
error_diff = zeros(new_frames,2);
for frame = 1:new_frames

    % Generate data
    [TX_1,TX_2,x_DD] = gen_data(bit_order,S,syms_per_f);
    TX_bit = Gamma_MN' * TX_1;
    TX_sym = Gamma_MN' * TX_2;
    x_tilde = Gamma_MN' * x_DD;

    % Generate channel
    t_offset = max_timing_offset * Ts;
    [chn_g,chn_tau,chn_v] = channel_generation(Fc,vel);
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
    % H_tilde = gen_H_tilde(H,M,N,Lp,Ln,F_N);

    % %% MUSIC PART
    % Settings
    noiseless = false;
    frame_samps = 1;
    sample_start = 1;
    sample_rate = 1;
    eig_thresh = 3e-14;
    v_res = 0.01;

    % Unwrap matrix to use all NM possible elements per diagonal
    H_left = [H(1:Lp,N*M-Lp+1:N*M); zeros(N*M-Lp,Lp)];
    H_right = [zeros(N*M+Ln,-Ln); H(N*M+Ln+1:N*M,1:-Ln)];
    H_center = H;
    H_center(1:Lp,N*M-Lp+1:N*M) = 0;
    H_center(N*M+Ln+1:N*M,1:-Ln) = 0;
    H_unwrapped = [H_left H_center H_right];

    % Sample main diagonal and take every Nth sample
    % size_samp = min(N*M+Ln+1,N*M-Lp+1);
    size_samp = N*M;
    x = zeros(ceil(size_samp / sample_rate),Lp-Ln-1);
    count = 0;
    x_cell = cell(frame_samps,1);
    for j = 1:frame_samps
        if ~noiseless
            z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
        end
        for i = 2:Lp-Ln
            count = count + 1;

            if noiseless
                h_diag = diag(H_unwrapped,i-1);
            else
                h_diag = diag(H_unwrapped + z_tilde,i-1);
            end
            x_samp = h_diag(sample_start:sample_rate:size_samp);
            x(:,i-1) = x_samp;
        end
        x_cell{j} = x;
    end
    x = vertcat(x_cell{:});

    % Generate covariance and select noise subspace
    R_x = x * x';
    [U,D] = eig(R_x);

    % Estimate number of paths
    if noiseless
        p_est = max(rank(D),sum(D > eig_thresh,"all"));
    else
        p_est = rank(D) - 1;
    end

    % Perform eigenvalue decomposition to find noise space
    eig_vals = diag(D);
    [~, ind] = sort(eig_vals, 'descend');
    Us = U(:,ind);
    Ds = D(ind,ind);
    U_N = Us(:,(p_est+1):end);

    % Sweep through all possible Doppler values and generate costs
    v_max = 100 * ceil((vel * (1000/3600)*Fc) / physconst('LightSpeed') / 100);
    v_range = -v_max:v_res:v_max;
    power_vec = ((sample_start-1):sample_rate:(size_samp-1))';
    power_vec = repmat(power_vec,frame_samps,1);
    e_temp = exp(1j*2*pi*Ts) .^ power_vec;
    e = e_temp.^v_range;
    e = num2cell(e,1);
    d_sqrd = cellfun(@(a) norm(U_N' * a)^2,e);
    P_hat = 1 ./ d_sqrd;

    % Estimate doppler shifts
    [~,v_idx] = findpeaks(P_hat);
    if length(v_idx) ~= p_est
        % figure(1)
        % % Plot cost function (-P_hat)
        % plot(v_range,-P_hat)
        % xlabel("Estimated Doppler")
        % ylabel("Cost function")

        [~,v_idx] = maxk(P_hat,p_est);
    end
    if length(v_idx) ~= length(chn_v)
        v_idx = v_idx * ones(1,length(chn_v));
    end
    v_est = v_range(v_idx);

    % Estimate first Doppler value
    v_est_sorted = sort(v_est);
    chn_v_sorted = sort(chn_v);
    for i = 1:length(chn_v)
        try
            fprintf("(Frame %d/%d) Estimated v_%d = %f, true v_%d = %f\n",frame,new_frames,i-1,v_est_sorted(i),i-1,chn_v_sorted(i));
        catch
            1;
        end
    end
    % R_x_reconst = U * D * U';

    error_diff(frame,:) = abs(v_est_sorted - chn_v_sorted);

    % % if nnz(error_diff(frame,:) > 100) > 0
    % if 1
    %     figure(1)
    %     % Plot cost function (-P_hat)
    %     plot(v_range,-P_hat)
    %     xlabel("Estimated Doppler")
    %     ylabel("Cost function")
    % 
    %     1;
    % end

    % if nnz(error_diff(frame,:) > 2) > 0
    %     figure(1)
    %     % Plot cost function (-P_hat)
    %     plot(v_range,-P_hat)
    %     xlabel("Estimated Doppler")
    %     ylabel("Cost function")
    % 
    %     1;
    % end

    % %% Code for simulating OTFS system
    % 
    % % Generate noise
    % z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));
    % 
    % % Create receive vector
    % y_tilde = H_tilde * x_tilde + z_tilde;
    % 
    % % Iterative Detector - JRW
    % switch receiver_name
    %     case "CMC-MMSE"
    %         [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_CMC_MMSE(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters,R);
    %     case "MMSE"
    %         [x_hat,iters_vec(frame),t_RXiter_vec(frame),t_RXfull_vec(frame)] = equalizer_MMSE(y_tilde,H_tilde,Es,N0);
    %     otherwise
    %         error("Unsupported receiver for the simulated system!")
    % end
    % 
    % % Hard detection for final x_hat
    % dist = abs(x_hat.' - S).^2;
    % [~,min_index] = min(dist);
    % RX_sym = min_index.' - 1;
    % 
    % % Convert final RX_sym to RX_bit
    % RX_bit = bit_order(RX_sym+1,:);
    % 
    % % Error calculation
    % bit_error_vec = TX_bit ~= RX_bit;
    % sym_error_vec = TX_sym ~= RX_sym;
    % bit_errors(frame) = sum(bit_error_vec,"all");
    % sym_errors(frame) = sum(sym_error_vec,"all");
    % if sum(bit_error_vec,"all") > 0
    %     frm_errors(frame) = 1;
    % else
    %     frm_errors(frame) = 0;
    % end

end

% % Get parameters for throughput
% frame_duration = N * T;
% bandwidth_hz = M / T;
% 
% % Calculate BER, SER and FER
% metrics.BER = sum(bit_errors,"all") / (new_frames*syms_per_f*log2(M_ary));
% metrics.SER = sum(sym_errors,"all") / (new_frames*syms_per_f);
% metrics.FER = sum(frm_errors,"all") / (new_frames);
% metrics.Thr = (log2(M_ary) * syms_per_f * (1 - metrics.FER)) / (frame_duration * bandwidth_hz);
% metrics.RX_iters = mean(iters_vec);
% metrics.t_RXiter = mean(t_RXiter_vec);
% metrics.t_RXfull = mean(t_RXfull_vec);

%

fprintf("Mean abs error per path: ")
for i = 1:size(error_diff,2)
    fprintf("%f", mean(error_diff(:,i)))
    if i ~= size(error_diff,2)
        fprintf(", ")
    else
        fprintf("\n")
    end
end
fprintf("Mode abs error per path: ")
for i = 1:size(error_diff,2)
    fprintf("%f", mode(error_diff(:,i)))
    if i ~= size(error_diff,2)
        fprintf(", ")
    else
        fprintf("\n")
    end
end
fprintf("Max abs error per path: ")
for i = 1:size(error_diff,2)
    fprintf("%f", max(error_diff(:,i)))
    if i ~= size(error_diff,2)
        fprintf(", ")
    else
        fprintf("\n")
    end
end

figure(10)
for i = 1:size(error_diff,2)
    subplot(size(error_diff,2),1,i)
    histogram(error_diff(:,i))
    title(sprintf("Error distribution for ray %d",i))
end