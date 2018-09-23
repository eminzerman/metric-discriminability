function clsStruct = metricDiscriminability(objResults, subjResMOS, subjStd, noOfStim, noOfSubj)
% METRICDISCRIMINABILITY Checks the objective quality metric results in
%   terms of metric discriminability.
%
%   Check the following paper for details:
%     E. Zerman, G. Valenzise, and F. Dufaux, "An Extensive Performance 
%     Evaluation of Full-Reference HDR Image Quality Metrics", Quality 
%     and User Experience, volume 2, April 2017.
%
% ---------------------
%  Copyright (C) 2018 - Emin Zerman
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------

% --- Find statistical equivalence for subjective quality results ---
% Check input parameters and find 'stats'
if (size(subjResMOS,1) ~= 1) && (size(subjResMOS',1) ~= 1)
    subjResArr = subjResMOS;
    subjResMOS = mean(subjResArr,2);
    [~, sortIdx] = sort(subjResMOS);
    try
        [~, ~, stats]= anova1(double(subjResArr(sortIdx,:)'));
    catch
        stats = getSubjStats(subjResMOS, subjStd, size(subjResArr,1), size(subjResArr,2), sortIdx);
    end
elseif (size(subjResMOS,1) == 1) || (size(subjResMOS',1) == 1)
    % Subjective result is a vector, get other data
    [~, sortIdx] = sort(subjResMOS);
    stats = getSubjStats(subjResMOS, subjStd, noOfStim, noOfSubj, sortIdx);
end

% Find the ground truth
gt = logical(1 - multiple_comparison(stats, N));
Npos = sum(gt(:));
Ipos = gt == 1;
Nneg = N^2 - Npos;
Ineg = gt == 0;




% --- Find the classification rates ---
% Find statistical equivalence fo objective quality results
% Generate classification results for each selected objective metric
countIdx = 1;
for ind = 1:size(objResults,2)
    qualVec = objResults(sortIdx,ind);
    minVal = min(qualVec);
    maxVal = max(qualVec);
    
    difMat = repmat(qualVec, [1 length(qualVec)]) - repmat(qualVec', [length(qualVec) 1]);
    
    % Find ROC curve for each
    countIdx2 = 1;
    % Sweep for the metric difference range
    for ind2 = logspace(0, 3, 50)
        difTh = ((maxVal - minVal)/(ind2));
        difMap = abs(difMat) < difTh;

        % Assign to classess
        TP = not(gt) & not(difMap);
        FP = gt & not(difMap);
        FN = not(gt) & difMap;
        TN = gt & difMap;
        difMaps(:,:,countIdx2) = difMap;
        
        % Find the counts and classification rates
        tp = sum(TP(Ineg));
        fp = sum(FP(Ipos));
        tn = sum(TN(Ipos));
        fn = sum(FN(Ineg));
        tpr = tp/Nneg;
        fpr = fp/Npos;
        tnr = tn/Npos;
        fnr = fn/Nneg;

        % accuracy
        acc = (tp+tn)/(tp+tn+fp+fn);
        % balanced accuracy
        accB = ((0.5*tp)/(tp+fn))+((0.5*tn)/(tn+fp));
        % precision 
        prec = tp/(tp+fp);
        % recall == TPR == 
        recl = tp/(tp+fn);
        
        tpVec(countIdx2) = tp;
        fpVec(countIdx2) = fp;
        tnVec(countIdx2) = tn;
        fnVec(countIdx2) = fn;
               
        tprVec(countIdx2) = tpr;
        fprVec(countIdx2) = fpr;
        tnrVec(countIdx2) = tnr;
        fnrVec(countIdx2) = fnr;
        
        accVec(countIdx2) = acc;
        accBVec(countIdx2) = accB;
        precVec(countIdx2) = prec;
        reclVec(countIdx2) = recl;
        
        countIdx2 = countIdx2 + 1;
    end
    
    % Store ROC values
    ROC(1,:) = fprVec;
    ROC(2,:) = tprVec;
    accBal = accBVec;
    T = (maxVal - minVal)./logspace(0, 4, 70);
    fprVec2 = fprVec; fprVec2(end+1) = 1;
    tprVec2 = tprVec; tprVec2(end+1) = 1;
    AUC = abs(trapz(fprVec2, tprVec2));
    MeanAvP = abs(trapz(precVec, reclVec));
    
    % Find optimal point
    [~, id] = min(sqrt((-fprVec).^2 + (1-tprVec).^2));
    optTh = T(id);
    optX = fprVec(id);
    optY = tprVec(id);

    % Store all th data    
    rocs(:,:,countIdx) = ROC;
    accBs(:,countIdx) = accBal;
    ts(:,:,countIdx) = T;
    aucs(:,:,countIdx) = AUC;
    optths(countIdx) = optTh;
    optPts(:,:,countIdx) = [optX; optY];
    names.metric(countIdx) = dbStruct.names.metric(ind);
    labels.metric(countIdx) = dbStruct.labels.metric(ind);
    
    [~, id] = min((fprVec-0.05).^2);
    tP005 = T(id);
    accP005 = accBal(id);
    [~, id] = min((fprVec-0.25).^2);
    tP025 = T(id);
    accP025 = accBal(id);
    [~, id] = max(accBal);
    tMaxAcc = T(id);
    accMaxAcc = accBal(id);
    result(:,countIdx) = [AUC; MeanAvP; tP005; accP005; tP025; accP025; tMaxAcc; accMaxAcc];
    
    countIdx = countIdx+1;
    
end

clsStruct = struct('T', ts,...
                   'AUC', aucs,...
                   'OptPt', optPts,...
                   'OptTh', optths,...
                   'names', names,...
                   'labels', labels,...
                   'result', result,...
                   'accBs', accBs,...
                   'ROCs', rocs);
end

% ===================
function diff_matrix = multiple_comparison(stats, N)

[cval, ~] = multcompare(stats);
for pair = 1:size(cval,1)
    if (cval(pair,3)<0 && cval(pair,5)<0)  || ...
            (cval(pair,3)>0 && cval(pair,5)>0 ) % Ci interval does not contain 0
        cval(pair,6) = 1; % mean are significantly different
    else
        cval(pair, 6) = 0; % mean are not significantly different
    end
end

multcomp_disp = zeros(N);
multcomp_disp(sub2ind([N N],cval(:,1), cval(:,2))) = cval(:,6);

% figure

diff_matrix = multcomp_disp + multcomp_disp';

% imagesc(diff_matrix), axis image, colormap gray
end

function stats = getSubjStats(subjMOS, subjStd, noOfStim, noOfSubj, sortIdx)
% Generate stats struct for multcompare
% This operation is done for databases without raw subject data

stimL = noOfStim;
subjN = noOfSubj;
    
% Get gnames
stats1.gnames = [];
for ind = 1:stimL
    stats1.gnames = cat(1, stats1.gnames, sprintf('%3d',ind));
end

% Get n
stats1.n = repmat(subjN, 1, stimL);
stats1.source = 'anova1';

% Get means
stats1.means = subjMOS(sortIdx)';

% Get df
stats1.df = (subjN - 1)*stimL;

% Get s
varVec = subjStd(sortIdx).^2;
sqErr = varVec*(subjN-1);
sse = sum(sqErr);
stats1.s =sqrt(sse/stats1.df);

% Return
stats = stats1;

end