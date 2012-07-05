% parameters
L = 600
dx = L./[30 60 120 180 240]
maxspeed = 181.366

% errors
maxvec = [0.6426 0.1610 0.0402 0.0179 0.0101]
avvec  = [0.40294 0.10071 0.02518 0.01119 0.00629]
prcntavvec = [0.22217 0.05553 0.01388 0.00617 0.00347]

loglog(dx,prcntavvec/100,'o')
axis([2 25 1e-5 1e-2])
set(gca,'XTick',[2.5 3.3333 5 10 20])
xlabel('\Delta x'),  ylabel('(average error)/(maximum speed)')


% >> testJerr
% L =
%    600
% dx =
%            20           10            5       3.3333          2.5
% maxspeed =
%        181.37
% maxvec =
%        0.6426        0.161       0.0402       0.0179       0.0101
% avvec =
%       0.40294      0.10071      0.02518      0.01119      0.00629
% prcntavvec =
%       0.22217      0.05553      0.01388      0.00617      0.00347
% >> polyfit(log(dx),log(prcntavvec),1)
% ans =
%        2.0002      -7.4963