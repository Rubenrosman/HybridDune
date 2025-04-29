% colors plot Coastal dynamcis
% see for example plot S1 V9.png

lijn_kleuren = colorcet('R1','N',15);          % version for 3 storms, 3 tides each
n_all = [1 2 3   7 8 9   13 14 15]
lijn_kleuren = lijn_kleuren(n_all,:);

%results:
lijn_kleuren = [
         0    0.1866    0.9626
    0.0948    0.3882    0.7161
    0.1883    0.5033    0.4919
    0.7730    0.7720    0.1165
    0.9433    0.7941    0.1423
    0.9790    0.7092    0.1287
    0.9769    0.3745    0.3615
    1.0000    0.4710    0.6640
    0.9923    0.5716    0.9805]
