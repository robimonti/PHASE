function mask=polynomialMask(count,originalDegree,protectedDegree,dimension)
%POLYNOMIALMASK Extra temporal terms must not introduce extra spatial terms.
mask=true(count,1);
if protectedDegree==0,return;end
if strcmp(dimension,'1D'),sizes=[1 3 6 10];timeIndices=[3 5];
else,sizes=[1 4 10 20];timeIndices=[4 7];end
mask(:)=false;mask(1:min(count,sizes(originalDegree+1)))=true;
mask(timeIndices(1:protectedDegree))=true;
end
