function [response_mat] = signrank_return_responsive_csbaseline(inputmat,pick_frame,frame_post,baseline_stop_frame,alpha)
%   Responsive: All cells with stimulus-evoked responses that were significantly different from baseline 
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
    
    response_mat=zeros(size(inputmat,1),size(pick_frame,2));
    for i=1:size(inputmat,1)
        for j=1:size(pick_frame,2)
            % ttest
            %response_mat(i,j)=ttest2(baseline(i,:),inputmat(i,pick_frame(j):pick_frame(j)+frame_post-1),"Alpha",alpha);
            % ranksum
            [p,~,~] = signrank(baseline(i,:,j), inputmat(i, pick_frame(j):pick_frame(j) + frame_post - 1));
            response_mat(i,j) = p <= alpha;
        end
    end

end