function [EVAL] = evaluate_scores(ACTUAL,PREDICTED)
% This fucntion evaluates the performance of a classification model by 
% calculating the common performance measures: Precision, Recall, F-Measure
% Input: ACTUAL = Column matrix with actual class labels of the training
%                 examples
%        PREDICTED = Column matrix with predicted class labels by the
%                    classification model
% Output: EVAL = Row matrix with all the performance measures
%

% Get the confusion matrix
conf_mat = confusionmat(ACTUAL, PREDICTED);

% Get number of classes
k = size(conf_mat,1);

% Initialize vectors
precision = zeros(k, 1);
recall    = zeros(k, 1);
F_score   = zeros(k, 1);

%% Compute precision, recall and F_score for each class
for i = 1:k
    tp = conf_mat(i,i);             % True Positive
    fp = sum(conf_mat(:,i)) - tp;   % False Positive
    fn = sum(conf_mat(i,:)) - tp;   % False Negative

    if tp + fp ~= 0
        precision(i) = tp / (tp + fp);
    else
        precision(i) = 0.0;
    end

    if tp + fn ~= 0
        recall(i) = tp / (tp + fn);
    else
        recall(i) = 0.0;
    end

    if precision(i) + recall(i) ~= 0
        F_score(i) = 2 * (precision(i) * recall(i)) /...
                         (precision(i) + recall(i));
    else
        F_score(i) = 0;
    end
end

%% Compute mean values
precision_avg   = mean(precision);
recall_avg      = mean(recall);
F_score_avg     = mean(F_score);

EVAL = [precision_avg recall_avg F_score_avg];

end