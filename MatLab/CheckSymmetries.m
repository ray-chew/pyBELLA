% check symmetry of field th() as generated while running Show_sequence.m
thsymm  = 0.5 * (th + th(:,end:-1:1));
thasymm = 0.5 * (th - th(:,end:-1:1));

figure(100)
subplot(2,1,1)
contour(th,)

