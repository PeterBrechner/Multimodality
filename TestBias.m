function TestBias(nSDs)

moments = [0 2];
binEdges = [linspace(0.005,0.12,116), linspace(0.13,1.29,117)];
bins = (binEdges(2:end)+binEdges(1:end-1))/2;
bins_diff = (binEdges(2:end)-binEdges(1:end-1));
alpha = 0.999;
total = 1e5;
sd = zeros(9*nSDs, length(bins));
sqError = zeros(9*nSDs, length(bins));
modality = strings(9*nSDs, 1);
trimodals = modality;
bimodal1s = modality;
bimodal2s = modality;
unimodals = modality;
N0s = zeros(9*nSDs,3);
mus = zeros(9*nSDs,3);
lambdas = zeros(9*nSDs,3);
meds = zeros(9*nSDs,3);
for ct = 1:4*nSDs
	trimodals(ct) = "Trimodal";
	bimodal1s(ct) = "Bimodal1";
	bimodal2s(ct) = "Bimodal2";
	unimodals(ct) = "Unimodal";
	mu = [-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5))];
	lambda = [108*(4+min(mu))*exp(normrnd(0,0.5)),...
		18*(4+median(mu))*exp(normrnd(0,0.5)),...
		6*(4+max(mu))*exp(normrnd(0,0.5))];
	ratios = exp(normrnd(-1,2,1,2));
	[sd(ct,:), sqError(ct,:), modality(ct), N0s(ct,:), mus(ct,:),...
		lambdas(ct,:), meds(ct,:)] =...
		BiasTest(moments, ratios, mu, lambda, alpha, binEdges, total);
end
for ct=(1+4*nSDs):6*nSDs
        mu = [-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5))];
	mu = sort(mu);
	mu(3) = mu(2);
        lambda = [108*(4+min(mu))*exp(normrnd(0,0.5)),...
                18*(4+median(mu))*exp(normrnd(0,0.5)),...
                6*(4+max(mu))*exp(normrnd(0,0.5))];
	lambda = flip(sort(lambda));
	lambda(3) = lambda(2);
	ratios = exp(normrnd(-1,2,1,2));
        [sd(ct,:), sqError(ct,:), modality(ct), N0s(ct,:), mus(ct,:),...
                lambdas(ct,:), meds(ct,:)] =...
                BiasTest(moments, ratios, mu, lambda, alpha, binEdges, total);
end
for ct=(1+6*nSDs):8*nSDs
        mu = [-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5))];
	mu = sort(mu);
        mu(1) = mu(2);
        lambda = [108*(4+min(mu))*exp(normrnd(0,0.5)),...
                18*(4+median(mu))*exp(normrnd(0,0.5)),...
                6*(4+max(mu))*exp(normrnd(0,0.5))];
	lambda = flip(sort(lambda));
        lambda(1) = lambda(2);
	ratios = exp(normrnd(-1,2,1,2));
        [sd(ct,:), sqError(ct,:), modality(ct), N0s(ct,:), mus(ct,:),...
                lambdas(ct,:), meds(ct,:)] =...
                BiasTest(moments, ratios, mu, lambda, alpha, binEdges, total);
end
for ct=(1+8*nSDs):9*nSDs
        mu = [-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5)),-1+3*exp(normrnd(0,0.5))];
        mu(3) = mu(2);
	mu(1) = mu(2);
        lambda = [108*(4+min(mu))*exp(normrnd(0,0.5)),...
                18*(4+median(mu))*exp(normrnd(0,0.5)),...
                6*(4+max(mu))*exp(normrnd(0,0.5))];
        lambda(3) = lambda(2);
	lambda(1) = lambda(2);
	ratios = exp(normrnd(-1,2,1,2));
        [sd(ct,:), sqError(ct,:), modality(ct), N0s(ct,:), mus(ct,:),...
                lambdas(ct,:), meds(ct,:)] =...
                BiasTest(moments, ratios, mu, lambda, alpha, binEdges, total);
end
time = linspace(1,9*nSDs,9*nSDs).';
fit_moments = [0 2 4];
iwc = ones(9*nSDs,1);
rawcount = sqError;
%rawcount(1:6) = rawcount(1:6).*bins(1:6).^2./bins(7);
%sqError(1:6) = sqError(1:6)*bins(7)/bins(1:6).^2;
%rawcount(76:end) = rawcount(76:end)*100;
%sqError(76:end) = sqError(76:end)/100;
sd = 1e-5*sd;
sqError = 1e-10*sqError;
[tr, ba, bb, un, hix, hix2, mp3, mp4, mp5, mp6] =...
    Input(time, fit_moments, bins, bins_diff, sd, iwc, rawcount, sqError);
mTr = find(strcmp(modality, "trimodal"));
mBa = find(strcmp(modality, "bimodal1"));
mBb = find(strcmp(modality, "bimodal2"));
mUn = find(strcmp(modality, "unimodal"));
fpa = intersect([tr; ba], [mBb; mUn]);
fpb = intersect([tr; bb], [mBa; mUn]);
fsa = intersect([ba], [mBb]);
fsb = intersect([bb], [mBa]);
fna = intersect([bb; un], [mTr; mBa]);
fnb = intersect([ba; un], [mTr; mBb]);
length(mTr)
length(mBa)
length(mBb)
length(mUn)
length(tr)
length(ba)
length(bb)
length(un)
length(fpa)
length(fpb)
length(fsa)
length(fsb)
length(fna)
length(fnb)
[mus(fpa,1), lambdas(fpa,1), mus(fpa,2), lambdas(fpa,2), mus(fpa,3),...
 lambdas(fpa,3)];
[mus(fpb,1), lambdas(fpb,1), mus(fpb,2), lambdas(fpb,2), mus(fpb,3),...
 lambdas(fpb,3)];
[mus(fna,1), lambdas(fna,1), mus(fna,2), lambdas(fna,2), mus(fna,3),...
 lambdas(fna,3)];
[mus(fnb,1), lambdas(fnb,1), mus(fnb,2), lambdas(fnb,2), mus(fnb,3),...
 lambdas(fnb,3)];
%[un, mp6(:,2), mus(un,2), mp6(:,3), lambdas(un,2)]
%modality
PlotMedians2(1e4*bins, 1e3*sd(mTr,:), 1e3*sd(mBa,:), 1e3*sd(mBb,:), 1e3*sd(mUn,:), mTr, mBa, mBb,...
 mUn, iwc, "g/m^3", 'HIWCimagesLatest/actual');
PlotMedians2(1e4*bins, 1e3*sd(fpb,:), 1e3*sd(fnb,:), 1e3*sd(fpa,:), 1e3*sd(fna,:), fpb, fnb, fpa,...
 fna, iwc, "g/m^3", 'HIWCimagesLatest/false');

figure
scatter(meds(un,1)./meds(un,2),meds(un,2)./meds(un,3),25,'b','filled')
hold on
scatter(meds(ba,1)./meds(ba,2),meds(ba,2)./meds(ba,3),25,'m','filled')
scatter(meds(bb,1)./meds(bb,2),meds(bb,2)./meds(bb,3),25,'g','filled')
scatter(meds(tr,1)./meds(tr,2),meds(tr,2)./meds(tr,3),25,'y','filled')
hold off
xlabel('Mid Over Small Median Diameter')
ylabel('Large Over Mid Median Diameter')
title('Median Diameter Ratios')
saveas(gca,'HIWCimagesLatest/MDRatios.png')
