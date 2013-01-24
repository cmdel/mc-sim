% Generate Quasi-random
% Uniform variates
n = 1000;
genid                            = int64(1);
idim                             = int64(2);  % On Uv for the QEvariance
iskip                            = int64(1000);

rcord                            = int64(1);

[iref,  ifail]             = g05yl(genid, idim, iskip);

[quasi, iref, ifail]             = g05ym(int64(n), rcord, iref); 

x = rand(1,1000);
y = rand(1,1000);

box on;
figure(1)
subplot 121
scatter(quasi(1,:),quasi(2,:));
subplot 122
scatter(x,y);

[iref, ifail] = g05yl(int64(3), idim, iskip);
[quasi, iref, ifail] = g05ym(int64(n), rcord, iref);
figure(2)
subplot 121
scatter(quasi(1,:),quasi(2,:));
subplot 122
scatter(x,y);

[iref, ifail] = g05yl(int64(4), idim, iskip);
[quasi, iref, ifail] = g05ym(int64(n), rcord, iref);
figure(3)
subplot 121
scatter(quasi(1,:),quasi(2,:));
subplot 122
scatter(x,y);
  
