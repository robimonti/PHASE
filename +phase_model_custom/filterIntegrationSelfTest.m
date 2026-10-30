function result=filterIntegrationSelfTest()
%FILTERINTEGRATIONSELFTEST Exercise every custom entry through geoSplinter.
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
c.step_t_ST=30;
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
        [~,outDates,~,ps,model,stds]=feval(['STmodel_' family{1} dimension{1} '_custom'],args{:},c,folder);
        assert(numel(outDates)==numel(dates) && size(ps,1)==4 && size(ps,2)==numel(dates));
        assert(all(isfinite(model(:))) && all(isfinite(stds(:))) && all(model(:,1)==0));
        diagnostic=load(fullfile(folder,'custom_model_results.mat'),'clean');
        if isempty(firstClean), firstClean=diagnostic.clean;
        else, assert(isequaln(firstClean,diagnostic.clean),'Outlier filtering differs across custom entries.'); end
    end
end
result=struct('ok',true,'message','All four custom entries passed the native outlier-filter integration test.');
disp(result.message);
end
