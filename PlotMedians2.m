function PlotMedians1(bins, trimodalSds, bimodal1Sds, bimodal2Sds,...
    unimodalSds, trimodalIndices, bimodal1Indices, bimodal2Indices,...
    unimodalIndices, normalizeBy, normUnits, dir)

%bins = row vector containing each bin's center diameter, in um
%trimodalSds = subset of trimodal SDs, for example with IWC > 1.5 g/m^3...
    %or with T < -40C, in #/L/um
%bimodal1Sds = subset of bimodal1 SDs, for example with IWC > 1.5 g/m^3...
    %or with T < -40C, in #/L/um
%bimodal2Sds = subset of bimodal2 SDs, for example with IWC > 1.5 g/m^3...
    %or with T < -40C, in #/L/um
%unimodalSds = subset of unimodal SDs, for example with IWC > 1.5 g/m^3...
    %or with T < -40C, in #/L/um
%trimodalIndices = indices of normalizeBy vector corresponding to trimodalSds
    %trimodalIndices = column vector
%bimodal1Indices = indices of normalizeBy vector corresponding to bimodal1Sds
    %bimodal1Indices = column vector
%bimodal2Indices = indices of normalizeBy vector corresponding to bimodal2Sds
    %bimodal2Indices = column vector
%unimodalIndices = indices of normalizeBy vector corresponding to unimodalSds
    %unimodalIndices = column vector
%normalizeBy = vector to normalize sds by, for example IWC, concentration,...
    %moment 2 of concentration, etc
%normUnits = string containing units of normalizeBy vector, for example...
    %'g/m^{3}' for normalizeBy = IWC

tri = trimodalIndices; %shorthand
b1i = bimodal1Indices; %shorthand
b2i = bimodal2Indices; %shorthand
uni = unimodalIndices; %shorthand


    if length(tri) ~= 0
        trnorm = trimodalSds./repmat(normalizeBy(tri),1,length(bins));
    else
        trnorm = trimodalSds;
    end
    if length(tri) ~= 1
	trnorm = median(trnorm);
    end
    if length(b1i) ~= 0
        b1norm = bimodal1Sds./repmat(normalizeBy(b1i),1,length(bins));
    else
        b1norm = bimodal1Sds;
    end
    if length(b1i) ~= 1
	b1norm = median(b1norm);
    end
    if length(b2i) ~= 0
        b2norm = bimodal2Sds./repmat(normalizeBy(b2i),1,length(bins));
    else
        b2norm = bimodal2Sds;
    end
    if length(b2i) ~= 1
	b2norm = median(b2norm);
    end
    if length(uni) ~= 0
        unnorm = unimodalSds./repmat(normalizeBy(uni),1,length(bins));
    else
        unnorm = unimodalSds;
    end
    if length(uni) ~= 1
	unnorm = median(unnorm);
    end


figure
plot(bins.', 1e9*unnorm.', 'Color', [0 0.7 1], 'DisplayName', 'Unimodal');
hold on
plot(bins.', 1e9*b1norm.', 'm', 'DisplayName', 'Bimodal1');
plot(bins.', 1e9*b2norm.', 'g', 'DisplayName', 'Bimodal2');
plot(bins.', 1e9*trnorm.', 'Color', [1 0.7 0], 'DisplayName', 'Trimodal');
hold off
set(gca,'xscale','log');
set(gca,'yscale','log');
axis([10 20000 1e3 1e13]);
xlabel('Diameter [{\mu}m]')
ylabel(strcat('Normalized Conc [m^{-4}/(', normUnits, ')]'))
title('Conc and Normalized Median Modality Curves')
legend
saveas(gca,strcat(dir,'Norm_Med_Modality_Curves.png'))
