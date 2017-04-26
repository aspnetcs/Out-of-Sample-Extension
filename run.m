% re-produces all examples from paper

% figure 2 and 4
example1( 'plot', 'MNIST', 10000, 2.^[3:8] , [], inf, [], []);
example1( 'plot', 'poker', 10000, 2.^[3:8] , [], inf, [], []);
 
% figure 3, e = table 4
[~] = example2(1000, 10, 2.^[-15:-5], [1e-8]);
[~] = example2(1000, 10, [1], 2.^[-15:-5]);
[e] = example2(1000, 10, [0.1, 0.01, 0.001, 0.0001]);
 
% table 5
[ time_ns,  std_ns, time_ms, std_ms] = example3( [1000, 2000, 4000, 8000, 10000], 1000, 5000, [50,100,200,400,800], 10 , {'10', '20', '21', '22'});
 
% err_[?] = table 6
[ err_yeast] = example4(10, 'yeast', 1400, [50] , 200, 100, {'no_update', 'nystrom', 'full', '20','21','22'},  {'pertrubation_correction'});
[ err_mnist] = example4(10, 'MNIST', 5000, [50] , 50, 100, {'no_update', 'nystrom', 'full', '20','21','22'},  {'pertrubation_correction'});
[ err_poker] = example4(10, 'poker', 5000, [50] , 50, 100, {'no_update', 'nystrom', 'full', '20','21','22'},  {'pertrubation_correction'});

 
