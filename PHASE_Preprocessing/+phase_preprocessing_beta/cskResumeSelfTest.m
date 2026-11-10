function cskResumeSelfTest()
%CSKRESUMESELFTEST Real provenance comparisons, including old 3-digit files.
folder=tempname;mkdir(fullfile(folder,'graphs'));cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
bounds=[9.77528 9.8105398 45.0586109 45.0775897];
old=[9.718 9.823 45.044 45.087];writeGraph('CSK_subsetMgraph2run.xml',old);writeGraph('CSK_subsetgraph2run.xml',old);
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeAreaChanged');
phase_preprocessing_beta.validateCskResume(folder,old,1,3);
% A full rebuild is allowed despite old provenance.
phase_preprocessing_beta.validateCskResume(folder,bounds,0,1);
% Unchanged AOIs from releases that rounded to 3 decimals remain resumable.
writeGraph('CSK_subsetMgraph2run.xml',round(bounds,3));writeGraph('CSK_subsetgraph2run.xml',round(bounds,3));
phase_preprocessing_beta.validateCskResume(folder,bounds,1,3);
writeGraph('CSK_subsetgraph2run.xml',old);
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeAreaChanged');
delete(fullfile(folder,'graphs','CSK_subsetMgraph2run.xml'));
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeProvenanceMissing');
fprintf('CSK resume provenance regressions passed.\n');
    function writeGraph(name,b)
        fid=fopen(fullfile(folder,'graphs',name),'w');
        fprintf(fid,'<graph><geoRegion>POLYGON ((%.10f %.10f,%.10f %.10f,%.10f %.10f,%.10f %.10f,%.10f %.10f))</geoRegion></graph>', ...
            b(1),b(3),b(2),b(3),b(2),b(4),b(1),b(4),b(1),b(3));fclose(fid);
    end
end
function mustFail(action,id)
try,action();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.message);return;end
error('Expected failure: %s',id);
end
