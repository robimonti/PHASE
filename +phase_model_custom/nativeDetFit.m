function residual=nativeDetFit(c,r,ids,fit,config)
%NATIVEDETFIT Original spline selection/fit; recover its actual linear map.
root=phase_model_beta.projectRoot(); addpath(fullfile(root,'MatlabFunctions'));
state=struct(); suffix=['DET' config.projDim];
state.varNoise_method=config.(['varNoise_' suffix]);
state.varNoise_manual=config.(['varNoise_manual_' suffix]);
state.num_spl_method=config.(['num_spl_method_' suffix]);
state.spline_method=config.(['spline_method_' suffix]);
state.num_spl_row_manual=config.(['num_spl_row_manual_' suffix]);
state.num_spl_col_manual=config.(['num_spl_col_manual_' suffix]);
state.num_spl_time_manual=NaN;
if strcmp(config.projDim,'2D'), state.num_spl_time_manual=config.num_spl_t_manual_DET2D; end
state.lambda_method=config.(['lambda_method_' suffix]);
state.lambda_manual=config.(['lambda_manual_' suffix]);
state.t_relIN=fit.context.times; state.xyIN_AOI=fit.context.xy;
state.nt=numel(state.t_relIN); state.np=size(state.xyIN_AOI,1);
state.alpha_sig=0.05;
[~,it]=ismember(c(:,4),state.t_relIN);
linear=sub2ind([state.np,state.nt],ids,it);
matrix=nan(state.np,state.nt); matrix(linear)=r;
state.data_final_clean=matrix;
state.res_spl=fit.context.filterResidual;
state.phi_3=@(v) (v>=0 & v<1).*((2-v).^3-4*(1-v).^3)/6+(v>=1 & v<2).*(2-v).^3/6;
state.phi_4=@(v)(v>=0 & v<2).*(2-v).^3/6;
base=fullfile(fit.context.outputDir,'native_residual');
for name={'input','output','jobs'}
    folder=fullfile(base,name{1}); if ~isfolder(folder), mkdir(folder); end
end
state.gS_input_path=fullfile(base,'input');state.gS_output_path=fullfile(base,'output');state.gS_job_path=fullfile(base,'jobs');
state.gS_dir=fullfile(root,'geoSplinter',detectOS());
state.s_SAR=fit.context.s;
state.grid_s_flt_col=reshape(repmat(state.s_SAR,1,state.nt)',[],1);
state.grid_t_flt_col=reshape(repmat(state.t_relIN,state.np,1)',[],1);
mask=~isnan(reshape(matrix',[],1));
state.grid_s_flt_col=state.grid_s_flt_col(mask);state.grid_t_flt_col=state.grid_t_flt_col(mask);
if strcmp(config.projDim,'1D'), native=phase_model_custom.nativeDet1D(state);
else, native=phase_model_custom.nativeDet2D(state); end
basis=str2func(['phase_model_custom.nativeBasis' config.projDim]);
B=basis(c,native);
if strcmp(config.projDim,'2D')
    N=B'*B+native.lambda_spl*speye(size(B,2));
    if condest(N)>1e10, coefficientMap=pinv(full(N))*B'; else, coefficientMap=N\B'; end
    coefficients=native.coef_spl;
else
    coefficients=native.nativeCoefficients(:);
    % geoSplinter uses its own Laplacian regularisation. Recover that exact
    % operator on the observation design's column space, rather than replace
    % it with a ridge penalty or invert the rounded diagnostic .mat.txt file.
    [U,S,~]=svd(full(B),'econ'); count=nnz(diag(S)>max(size(B))*eps(max(diag(S))));
    U=U(:,1:count); responses=zeros(size(B,2),count);
    pair=native.row_col_pairs(native.var_idx,:);
    [~,order]=sortrows([ids,c(:,4)],[1,2]);
    for k=1:count
        phase_model_custom.checkStop(config);
        writematrix([c(order,3),c(order,4),U(order,k)],fullfile(state.gS_input_path,'operator.txt'),'Delimiter','space');
        jobFile_analysis2D(2,2,'operator.txt','operator',numel(r),pair(1),pair(2), ...
            native.knots_x(2),native.knots_x(end-1),native.knots_y(end-1),native.knots_y(2), ...
            native.lambda_spl,14,state.gS_input_path,state.gS_output_path,state.gS_job_path);
        [status,output]=runGeoSplinter(fullfile(state.gS_dir,'geoSplinter_analysis'),fullfile(state.gS_job_path,'operator.job'));
        if status~=0, error('PHASE:customNativeOperator','geoSplinter operator failed: %s',output); end
        response=load(fullfile(state.gS_output_path,'operator.par.txt'));responses(:,k)=response(:);
    end
    coefficientMap=responses*U';
end
if norm(coefficientMap*r-coefficients)>1e-4*max(1,norm(coefficients))
    error('PHASE:customNativeOperator','Native spline sensitivity does not reproduce the native solution.');
end
residual=struct('kind','DET','native',true,'nativeState',native, ...
    'basis',basis,'coefficientMap',coefficientMap,'coefficients',coefficients, ...
    'noiseVariance',native.varNoise,'lambda',native.lambda_spl);
end
