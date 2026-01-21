function [x_hat,iters,t_RXiter,t_RXfull] = equalizer_SIC_MMSE(y,H)
% SIC-MMSE receiver seen in:
% "Iterative MMSE Detection for Orthogonal Time Frequency Space Modulation"
%     by Dr. Jinhong Yuan and Dr. Hai Lin
%
% Coded by JRW, 1/21/2026

% % Start runtime
% tStartRX = tic;

