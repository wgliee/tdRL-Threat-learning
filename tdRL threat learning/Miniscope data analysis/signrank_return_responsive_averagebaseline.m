function [response_mat] = signrank_return_responsive_averagebaseline(inputmat,pick_frame,frame_post,baseline_stop_frame,alpha)
%   Responsive: All cells with TRAIL-AVERAGED stimulus-evoked responses that were significantly different from baseline 
%   activity (significance criterion: P≤0.01) were classified as CS- or
%   US-responsive.
%   input: n*t
%       returns a response mat 
    if nargin<5
        alpha=0.01;
    end

    % create baseline for each stimulus shape:n*t*trail
    baseline = zeros(size(inputmat,1),frame_post,size(pick_frame,2));
    for i=1:size(pick_frame,2)
        baseline(:,:,i) = inputmat(:, baseline_stop_frame(i)-frame_post+1 : baseline_stop_frame(i));
    end
    % trail-averaged bl
    baselline_use = mean(baseline,3);
    
    response_raw = zeros(size(inputmat,1),frame_post,size(pick_frame,2));
    % n*t*trails
    for j=1:size(pick_frame,2)
        response_raw(:,:,j) = inputmat(:, pick_frame(j):pick_frame(j) + frame_post - 1);
    end
    response_use = mean(response_raw,3);

    for i=1:size(response_use,1)
            % ranksum
            [p,~,~] = signrank(baselline_use(i,:), response_use(i,:));
            response_mat(i) = p <= alpha;
    end

end