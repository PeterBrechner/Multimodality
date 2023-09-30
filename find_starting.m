function [starting] = find_starting(mag, coff1, coff2, Dmin, Dmax)

% for one SD
starting = 0;
mag

if coff1 == Dmin & coff2 == Dmax
    starting = [mag(1), 2, 60];
end

mu1 = min(0.4+mag(2), 0.6+0.2*mag(2));
la1 = 6/coff1;
mu2 = 1+mag(2);
la2 = find_starting_la(mag(3),mu2-mag(2),coff1,coff2);
mu3 = 4.2+0.6*mag(2);
la3 = 0.5*la2;
lower = find_starting_M1(mu1, coff1*la1, mu2, coff1*la2);
upper = find_starting_M1(mu3, coff2*la3, mu2, coff2*la2);

if coff1 == Dmin & coff2 ~= Dmax
    starting = [mag(1)/(1+upper), mu2, la2,...
        mag(1)*upper/(1+upper), mu3, la3];
end

if coff1 ~= Dmin & coff2 == Dmax
    starting = [mag(1)*lower/(1+lower), mu1, la1,...
        mag(1)/(1+lower), mu2, la2];
end

if coff1 ~= Dmin & coff2 ~= Dmax
    starting = [mag(1)*lower/(1+lower+upper), mu1, la1,...
        mag(1)/(1+lower+upper), mu2, la2,...
        mag(1)*upper/(1+lower+upper), mu3, la3];
end
