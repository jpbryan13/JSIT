function [qq_t,qqc_t,t,m,bfct,ctb,ctc] = getThresholdHist(qq,nA,nI,nD,nC,nB,tgt)
% getThresholdHist performs adaptive filtering, attempting to reduce the 
% number of false positives by removing transcript calls that appear
% similar to "blank barcodes" 
%Inputs:
%   qq: List of transcript calls
%   nA: Number of bins in cluster size histogram
%   nI: Number of bins in image intensity histogram
%   nD: Number of bins in X matrix value histogram
%   nC: Number of coding barcodes in the codebook
%   nB: Number of blank barcodes in the codebook
%   tgt: Target misidentification rate
%Outputs:
%   qqt: Filtered list of transcript calls, including blank barcodes
%   qqct: Filtered list of transcript calls, including only coding barcodes
%   t: Blank barcode fraction threshold to approach target mis-ID rate
%   m: mis-ID rate achieved by thresholding by t
%   bfct: Distribution of fraction of blank barcodes, in histogram of
%   cluster size, image intensity, and X matrix value
%   ctb: histogram of coding barcodes, according to above properties
%   ctc: histogram of blank barcodes, according to above properties

%Get data to assemble histograms 
aix = qq(:,end-2:end);
aix(:,2) = aix(:,2)-min(aix(:,2));
aix(:,3) = aix(:,3)-min(aix(:,3));
aix(:,2:3) = log10(aix(:,2:3)+1);

%Get histogram edges
aEdges = 1:nA;
iEdges = linspace(min(aix(:,2)),max(aix(:,2)),nI);
xEdges = linspace(min(aix(:,3)),max(aix(:,3)),nD);

%Create 3-dimensional histogram of all transcript calls 
[ct,~,~,loc] = histcn(aix,aEdges,iEdges,xEdges);
baix = aix(qq(:,end-3)>nC,:); %data for blank barcodes
caix = aix(qq(:,end-3)<=nC,:); %data for coding barcodes

if(~isempty(baix))    
    [ctb,~,~,locb] = histcn(baix,aEdges,iEdges,xEdges); %3D histogram for blanks
    [ctc,~,~,locc] = histcn(caix,aEdges,iEdges,xEdges); %3D histogram for coding barcodes
    %Correct for mismatched dimensions
    if(size(ctb,1)~=nA)
        ctb = padarray(ctb,[nA-size(ctb,1),0,0],'post');
    end
    if(size(ctc,1)~=nA)
        ctc = padarray(ctc,[nA-size(ctc,1),0,0],'post');
    end
    if(sum(size(ctb)~=size(ctc))~=0)
        try
            ctb = padarray(ctb,[size(ctc,1)-size(ctb,1),0,0],'post');
        catch
            ctc = padarray(ctc,[size(ctb,1)-size(ctc,1),0,0],'post');
        end
        try
            ctb = padarray(ctb,[0,size(ctc,2)-size(ctb,2),0],'post');
        catch
            ctc = padarray(ctc,[0,size(ctb,2)-size(ctc,2),0],'post');
        end
        try
            ctb = padarray(ctb,[0,0,size(ctc,3)-size(ctb,3)],'post');
        catch
            ctc = padarray(ctc,[0,0,size(ctb,3)-size(ctc,3)],'post');
        end
    end
    %Compute fraction of blank barcodes in each histogram bin
    bfct = ctb./(ctc+ctb);
    
    %Compute misidentification rate, filtering by various blank fraction
    %thresholds
    ts = 0.01:0.001:0.1;
    mRate = zeros(length(ts),1);
    for x = 1:length(ts)
        bst = sum(ctb(bfct<ts(x)));
        cst = sum(ctc(bfct<ts(x)));
        mRate(x) = (bst./nB)./(cst./nC);
    end
    
    %Find blank fraction threshold which gives mis-ID rate closest to
    %target
    [~,minInd] = min(abs(mRate-tgt));
    m = mRate(minInd);
    
    %Bookkeeping steps
    qqc = qq(qq(:,end-3)<=nC,:);
    inds_nz = (locc(:,1)~=0).*(locc(:,2)~=0).*(locc(:,3)~=0);
    locc_nz = locc(find(inds_nz),:);
    idx_nz = sub2ind([nA,nI,nD],locc_nz(:,1),locc_nz(:,2),locc_nz(:,3));
    
    %Filter list of coding transcript calls to achieve target mis-ID rate
    idx_t = find(bfct<ts(minInd));
    nzt = ismember(idx_nz,idx_t);
    finds_nz = find(inds_nz);
    finds_nzt = finds_nz(nzt);
    qqc_t = qqc(finds_nzt,:);
    
    %Filter list of all transcript calls to achieve target mis-ID rate
    inds_nz = (loc(:,1)~=0).*(loc(:,2)~=0).*(loc(:,3)~=0);
    loc_nz = loc(find(inds_nz),:);
    idx_nz = sub2ind([nA,nI,nD],loc_nz(:,1),loc_nz(:,2),loc_nz(:,3));
    idx_t = find(bfct<ts(minInd));
    nzt = ismember(idx_nz,idx_t);
    finds_nz = find(inds_nz);
    finds_nzt = finds_nz(nzt);
    qq_t = qq(finds_nzt,:);

    t = ts(minInd);
    else
        qq_t = qq;
        qqc_t = qq;
        t = 0;
        m = 0;
    end
end

