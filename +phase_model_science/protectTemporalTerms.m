function statistics=protectTemporalTerms(statistics,mask,degree,dimension)
%PROTECTTEMPORALTERMS Retain explicitly requested time terms during reduction.
% Indices match the original STmodel polynomial ordering; no new basis.
if degree==0,return;end
if strcmp(dimension,'1D'),indices=[3 5];else,indices=[4 7];end
protected=indices(1:degree);
statistics(ismember(find(mask),protected))=Inf;
end
