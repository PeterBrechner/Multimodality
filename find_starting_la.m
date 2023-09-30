function [la] = find_starting_la(lambda, mu_adj, coff1, coff2)

la = lambda+mu_adj*log(coff2/coff1)/(coff2-coff1);
