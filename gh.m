    X_pt=-2.88:0.01:2.88;
    Y_pt=(X_pt.^2)/4/f;
    plot(X_pt,Y_pt,'b')
    hold on
    alpha=0:pi/20:2*pi;
    X_abs=0.035*cos(alpha);
    Y_abs=0.035*sin(alpha)+1.71;
    plot(X_abs,Y_abs,'r')
    hold on