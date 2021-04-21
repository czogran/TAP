
G=transmit;
K=rss(2,2,2);
W=rss(2,6,2);

K=[pid(1) 1;
    1 pid(1)];



feedin = [1 2];
feedout = [1 2];
sys = feedback(G,K,feedin,feedout,-1);