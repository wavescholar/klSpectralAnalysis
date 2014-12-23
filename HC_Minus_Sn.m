% How much space is in a hypercube minus the unit ball?
% n= 1:1:20;
% D = HC_Minus_Sn(n);
% plot((D));

function dA = HC_Minus_Sn(n)

dA = 2.^ n -( pi .^ (n/2)/gamma( n/2 + 1) );
