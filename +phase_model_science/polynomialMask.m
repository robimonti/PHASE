function mask=polynomialMask(count,originalDegree,protectedDegree,dimension,design)
%POLYNOMIALMASK Retain native temporal interactions without extra spatial terms.
mask=true(count,1);
if protectedDegree==0,return;end
if strcmp(dimension,'1D')
    sizes=[1 3 6 10];linear=[3 6 9];quadratic=[5 10];
else
    sizes=[1 4 10 20];linear=[4 9 10 15 17 20];quadratic=[7 18 19];
end
original=min(count,sizes(originalDegree+1));
mask(:)=false;mask(1:original)=true;
extra=linear;if protectedDegree==2,extra=[extra quadratic];end
extra=sort(extra(extra<=count & extra>original));
% New interactions can be dependent for sparse/degenerate PS geometries.
% Keep the original basis untouched, admitting only independent additions.
for index=extra
    if nargin>=5
        A=design(:,[find(mask)' index]);
        A=A./max(vecnorm(A),eps);
        if rank(A,1e-9)<size(A,2),continue;end
    end
    mask(index)=true;
end
end
