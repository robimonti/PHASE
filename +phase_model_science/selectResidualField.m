function [state,B,M]=selectResidualField(c,r,ids,dimension,referenceDay,config)
%SELECTRESIDUALFIELD Select a low-dimensional residual mean, not a covariance.
% Held-out contiguous time blocks AND spatial blocks guard against absorbing
% stochastic fluctuations. Selection is conditional on the upstream custom
% fit. Fold scores are descriptive prediction diagnostics, not iid tests.
% M is the exact linear map for the selected weighted, regularised mean fit.
if nargin<6,config=struct();end
addpath(fullfile(phase_model_beta.projectRoot(),'MatlabFunctions'));
if strcmp(dimension,'1D'),z=c(:,3);else,z=c(:,1:2);end
[point,first,group]=unique(ids);positions=z(first,:);
origin=mean(positions,1);spatialScale=max(range(positions,1),1);
state=struct('kind','space-time residual field','dimension',dimension,'degree',0, ...
    'spatialDegree',0,'spatialKind','polynomial','spatialIntervals',0, ...
    'spatialBounds',[min((positions-origin)./spatialScale);max((positions-origin)./spatialScale)],'origin',origin,'spatialScale',spatialScale, ...
    'timeScale',max(range(c(:,4)),1),'coefficients',[],'lambda',0, ...
    'policy','Automatic residual field: blocked temporal and spatial validation', ...
    'selectionScope','Conditional on the upstream custom fit and retained outlier mask; not proof of stationarity.');
B=zeros(numel(r),0);M=zeros(0,numel(r));
[dates,~,epoch]=unique(c(:,4));
if numel(dates)<10||numel(point)<6
    state.reason='Insufficient independent locations/epochs for blocked field selection.';return
end
% Spatial folds are contiguous along the first principal axis, not random PS.
[~,~,V]=svd((positions-origin)./spatialScale,'econ');
projection=((positions-origin)./spatialScale)*V(:,1);
[distinctProjection,~,projectionGroup]=unique(projection);
spatial=min(5,ceil(projectionGroup*5/numel(distinctProjection)));
distinctPositions=size(unique(positions,'rows'),1);
folds=[min(5,ceil(epoch*5/numel(dates))),spatial(group)];
counts=accumarray(group,1);weights=1./counts(group);weights=weights/mean(weights);
base=zeros(2,5);
for axis=1:2
    for fold=1:5
        test=folds(:,axis)==fold;base(axis,fold)=weightedMSE(r(test),weights(test));
    end
end
records=struct('timeDegree',{},'spatialDegree',{},'spatialKind',{},'spatialIntervals',{},'lambda',{},'timeRatio',{},'spaceRatio',{},'score',{},'coefficientCount',{});
for timeDegree=1:2
    for candidate=1:6
        phase_model_custom.checkStop(config);
        phase_model_custom.reportProgress(config,54+((timeDegree-1)*6+candidate)/12, ...
            sprintf('Validating residual mean field: time degree %d, spatial candidate %d / 6',timeDegree,candidate));
        trial=state;trial.degree=timeDegree;spaceDegree=min(candidate-1,2);
        trial.spatialDegree=spaceDegree;
        if candidate>3
            if any(diff(state.spatialBounds,1,1)<=1e-10),continue;end
            trial.spatialKind='native cubic splines';trial.spatialIntervals=candidate-2;
        end
        A=phase_model_science.spaceTimeBasis(c,trial,referenceDay);p=size(A,2);
        % Spatial geometry and sample size must identify every retained term.
        if distinctPositions<max(6,2*p/timeDegree)||rank(A)<p,continue;end
        for lambda=[0 .001 .01 .1 1]
            phase_model_custom.checkStop(config);
            losses=nan(2,5);
            for axis=1:2
                for fold=1:5
                    test=folds(:,axis)==fold;train=~test;
                    if nnz(train)<4*p||rank(A(train,:))<p,continue;end
                    beta=solve(A(train,:),r(train),weights(train),lambda);
                    losses(axis,fold)=weightedMSE(r(test)-A(test,:)*beta,weights(test));
                end
            end
            if any(~isfinite(losses(:))),continue;end
            ratio=sum(losses,2)./max(sum(base,2),eps);
            records(end+1)=struct('timeDegree',timeDegree,'spatialDegree',spaceDegree, ...
                'spatialKind',trial.spatialKind,'spatialIntervals',trial.spatialIntervals, ...
                'lambda',lambda,'timeRatio',ratio(1),'spaceRatio',ratio(2), ...
                'score',max(ratio),'coefficientCount',p); %#ok<AGROW>
        end
    end
end
state.candidates=records;
if isempty(records),state.reason='No identifiable field passed blocked validation geometry.';return;end
acceptable=find([records.timeRatio]<.9 & [records.spaceRatio]<.9);
if isempty(acceptable),state.reason='No field improved both blocked errors by at least 10 percent.';return;end
best=min([records(acceptable).score]);
shortlist=acceptable([records(acceptable).score]<=1.1*best+1e-8);
% Prefer fewer coefficients; then lower held-out error. No user tuning needed.
[~,ordering]=sortrows([[records(shortlist).coefficientCount]' [records(shortlist).score]'],[1 2]);
chosen=records(shortlist(ordering(1)));
state.degree=chosen.timeDegree;state.spatialDegree=chosen.spatialDegree;state.spatialKind=chosen.spatialKind;
state.spatialIntervals=chosen.spatialIntervals;state.lambda=chosen.lambda;
B=phase_model_science.spaceTimeBasis(c,state,referenceDay);
[beta,M]=solve(B,r,weights,state.lambda);
if any(~isfinite(beta))
    state.degree=0;B=zeros(numel(r),0);M=zeros(0,numel(r));
    state.reason='Selected field was numerically ill-conditioned.';return
end
if sqrt(weightedMSE(B*beta,weights))<max(.5,.1*sqrt(weightedMSE(r,weights)))
    state.degree=0;B=zeros(numel(r),0);M=zeros(0,numel(r));
    state.reason='Validated correction was too small to justify additional mean terms.';return
end
state.coefficients=beta;state.selected=chosen;
state.reason='Selected the simplest field within 10 percent of the best joint blocked score.';
end
function [beta,M]=solve(A,y,w,lambda)
% Penalise scaled spatial/time coefficients, while retaining exact sensitivity.
rootWeight=sqrt(w/sum(w));
[Q,R]=qr([rootWeight.*A;sqrt(lambda)*eye(size(A,2))],0);
if rcond(R)<1e-8
    beta=nan(size(A,2),1);M=nan(size(A,2),numel(y));return
end
M=R\(Q(1:size(A,1),:)'.*rootWeight');beta=M*y;
end
function value=weightedMSE(y,w)
value=sum(w.*y.^2)/max(sum(w),eps);
end
