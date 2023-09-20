close all
clear all
import casadi.*

problem_options = NosnocProblemOptions();
problem_options.n_s = 2;
problem_options.dcs_mode = 'Step';
problem_options.print_level = 1;

model = NosnocModel();
x1 = SX.sym('x1');
model.x = x1;
model.c = [x1];
model.S = [-1; 1];

f_11 = 3;
f_12 = 1;
model.F = [f_11, f_12];
model.x0 = -1;

problem_options.T = pi/4;
problem_options.N_stages = 1;
problem_options.N_finite_elements = 2;
problem_options.cross_comp_mode = 7;
problem_options.lift_complementarities = 0;
model.dims.n_s = 2;

solver_options = NosnocSolverOptions();
solver_options.mpcc_mode = MpccMode.direct;
solver_options.lifting_phase1_tau = 0.0;

mpcc = NosnocMPCC(problem_options, model);
solver = NosnocSolver(mpcc, solver_options);
solver.set('alpha', [0,0,1,1]);
solver.set('lambda_p', [0,0,0.151,0.452]);
solver.set('lambda_n', [0.667,0,0,0]);
solver.set('x', [-0.667,0,0.151,0.452]);
solver.set('h', [0.3333,0.4521]);
if solver_options.mpcc_mode == MpccMode.lifting
    solver.set('c_lift', [nthroot(1.67,3),-nthroot(2,3),-nthroot(2,3),nthroot(0.603,3)]);
end

w0_base = solver.nlp.w0;
n_exp = 1000;
stds = [0.1,0.5,1];
converged = false(3*n_exp,1);
max_noise = zeros(3*n_exp,1);
max_infeasibilty = zeros(3*n_exp,1);
c=zeros(3*n_exp,3);
for j = 1:3
    noise_std = stds(j);
    for i = 1:n_exp
        k = i+(n_exp*(j-1));
        noise = noise_std*randn(size(w0_base));
        max_noise(k) = max(abs(noise));
        solver.set('w0', w0_base + noise);
        infeasibility = full(solver.nlp.g_fun(solver.nlp.w0, solver.p_val));
        max_infeasibility(k) = max(abs(infeasibility));
        [results,stats] = solver.solve();
        converged(k) = stats.converged();
        if converged(k)
            c(k,:) = [0,1,0];
        else
            c(k,:) = [1,0,0];
        end
    end
end
%%
figure('Position', [100, 100, 800, 700]);
scatter(max_noise, max_infeasibility,50,c, 'LineWidth', 1.2);
xlabel('$\ell_\infty$-norm of noise vector');
ylabel('$\ell_\infty$-norm of infeasiblity vector');
axis square
set(gca,'fontsize', 14);
if solver_options.mpcc_mode == MpccMode.lifting
    figure_name = ['lift_with_noise.pdf'];
    title (['Lifting'])
else
    figure_name = ['direct_with_noise.pdf']
    title (['Direct'])
end
exportgraphics(gca, figure_name, 'ContentType', 'vector');
