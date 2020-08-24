csvName = 'SettlingVelocities.csv';
opts = detectImportOptions(csvName);
T = readtable(csvName, opts);

names = T(:, 1);
logvelocities = T(:, 9);

figure(1)
scatter(logvelocities);
ylabel('log(velocity (m/s))')
xlabel('grain size')