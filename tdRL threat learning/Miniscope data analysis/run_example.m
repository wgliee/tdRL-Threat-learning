%% run_example.m
% Example script demonstrating how to run the main analyses included in
% the tdRL_analysis repository using the provided example dataset.
%
% The script loads "example.mat" and illustrates the expected workflow for:
%   1. Identification of neurons responsive to shock (US) and tone (CS)
%   2. Hierarchical clustering of US response
%   3. Functional connectivity / graph-related analysis
%   4. Population vector distance (PVD) analysis
%   5. Decoding analyses
%
% The example dataset is provided to demonstrate the required input data
% structure and basic usage of the analysis functions. 
%
% Required example data:
%   traceeach_zs  - z-scored neural activity traces organized as a cell array
%   freezing_each - freezing traces organized as a cell array
%   clustidx_each - neuronal cluster assignments organized as a cell array
%   traceeach_S   - deconvolved Ca2+ events used for high/low-freezing US decoding
%
% @Author:  Mingyang Wei, Shanghai Jiao Tong University, 2026
%% addpath
addpath("helper_functions\")
%% load data
clear
load("example.mat")
%% Responsiveness definition
pick_shock_frame=[210 360 480 660 750]*10;    %shock
shock_frame_post=40;
pick_tone_frame=[180 330 450 630 720]*10;    %tone
tone_frame_post=100;
baseline_stop_frame_tone = [180 330 450 630 720]*10 -1;  
baseline_stop_frame_shock = [210 360 480 660 750]*10 -1;  
alpha=0.01;

% responsiveness to stimuli
response_mat_shock = signrank_return_responsive_averagebaseline(traceeach_zs{1},pick_shock_frame,shock_frame_post,baseline_stop_frame_shock,alpha)';
response_mat_tone = signrank_return_responsive_averagebaseline(traceeach_zs{1},pick_tone_frame,tone_frame_post,baseline_stop_frame_tone,alpha)';
% responsiveness to stimuli for each trial
response_mat_shock_each = signrank_return_responsive_averagebaseline(traceeach_zs{1},pick_shock_frame,shock_frame_post,baseline_stop_frame_shock,alpha);
response_mat_tone_each = signrank_return_responsive_averagebaseline(traceeach_zs{1},pick_tone_frame,tone_frame_post,baseline_stop_frame_tone,alpha);
%% Hierarchial Clustering
tmp = traceeach_zs{1};
trace_shock_zs = [];
pre_baseline = 20;
frame_lag = 40;

for i = 1:length(pick_shock_frame)
    trace_shock_zs(:,:,i) = tmp(:,pick_shock_frame(i)+1-pre_baseline:pick_shock_frame(i)+frame_lag);
end
results = Hierarchial_Clustering(trace_shock_zs,4,pre_baseline,10);
%% graph related analysis
results = funcconnect_graph(traceeach_zs, freezing_each, clustidx_each);
%% PVD analysis
results_new = pvd_highlow_us_mahal(traceeach_zs, freezing_each, clustidx_each);
%% Decoding analysis
% decode US response
results = decoding_resp_function(traceeach_zs, "shock");

% decode high low freezing epoch
results_new = decode_high_low_freezing_cs_function(traceeach_zs,freezing_each,clustidx_each,'PlotResults',true);  %decode_high_low....

% train decoder on high/low freezing CS data, then apply the decoder to US data
results_new = hilow_decode_us([],[],traceeach_S,freezing_each,clustidx_each,'PlotResults',true);
