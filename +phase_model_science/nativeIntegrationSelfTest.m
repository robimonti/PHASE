function result=nativeIntegrationSelfTest()
%NATIVEINTEGRATIONSELFTEST Exercise optional controls in the four automatic engines.
root=phase_model_beta.projectRoot(); old=pwd; cd(root);
pathCleanup=onCleanup(@() cd(old)); %#ok<NASGU>
addpath(fullfile(root,'MatlabFunctions'));
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
input=fullfile(folder,'input'); output=fullfile(folder,'output');
jobs=fullfile(folder,'jobs'); synthesis=fullfile(folder,'synthesis');
mkdir(input); mkdir(output); mkdir(jobs); mkdir(synthesis);
xy=[500000 4980000;500040 4980000;500000 4980040;500040 4980040];
times=0:30:600; dates=datetime(2020,1,1)+days(times);
raw=repmat(0.03*(times-times(8)),4,1)+0.1*sin((1:4)'*(times+10)/40);
raw=raw-raw(:,8); ids=(1:4)';
center=struct('projected_distances',[0;40;0;40], ...
    'xy_centerline',xy(1:2,:),'s',[0;40]);
[xg,yg]=meshgrid([500000 500040],[4980000 4980040]);
aoi=xy([1,2,4,3,1],:);
c=phase_model_beta.defaultConfig(); c.deterministicModelMode='custom';
c.step_t_ST=30;c.temporalTrend='quadratic';
resultFolder=fullfile(folder,'files','customModel');
os=detectOS(); firstClean=[];
for dimension={'1D','2D'}
    c.projDim=dimension{1};
    for family={'DET','STC'}
        c.procType=['spatial' family{1}];
        switch [family{1} dimension{1}]
            case 'DET1D'
                args={raw,ids,dates,times,center,folder,6,input,output,jobs,synthesis, ...
                    os,'auto',NaN,'auto','MDL',NaN,NaN,'auto',NaN,'32 T',xy,30,'cleanObs',true,1,10,[]};
            case 'DET2D'
                args={raw,ids,dates,times,xg,yg,folder,6,input,output,jobs,os, ...
                    'auto',NaN,'auto','MDL',NaN,NaN,NaN,'auto',NaN,'32 T',xy,aoi,30,'cleanObs',true,1,10,[]};
            case 'STC1D'
                args={raw,ids,dates,times,center,folder,6,input,output,jobs,os, ...
                    NaN,NaN,'32 T',xy,30,'cleanObs',true,1,'gaussian','gaussian',10,[]};
            case 'STC2D'
                args={raw,ids,dates,times,xg,yg,true,folder,6,input,output,jobs,os, ...
                    NaN,NaN,'32 T',xy,aoi,30,'cleanObs',1,'gaussian','gaussian',10,[]};
        end
        args{end}=xy(1,:);
        figFolder=fullfile(folder,'figures');if ~isfolder(figFolder),mkdir(figFolder);end
        for k=1:numel(args),if ischar(args{k})&&strcmp(args{k},folder),args{k}=figFolder;end;end
        opts=phase_model_science.options(struct('temporalTrend','quadratic'));
        [~,outDates,~,ps,model,stds,query]=feval(['STmodel_' family{1} dimension{1}],args{:},opts);
        assert(numel(outDates)==numel(dates) && size(ps,1)==4);
        assert(all(isfinite(model(:)))&&all(isfinite(stds(:)))&&all(isfinite(query(:))));
        assert(isfolder(fullfile(folder,'files','residualDiagnostics')));
        fprintf('Automatic %s %s with residual controls passed.\n',family{1},dimension{1});
    end
end
result=struct('ok',true,'message','All four automatic entries completed with residual temporal controls.');
disp(result.message);
end
