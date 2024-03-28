function [results,stats,model,problem_options, solver_options] = test_simple_switch(irk_representation, irk_scheme, dcs_mode, cross_comp_mode)
    import casadi.*
    problem_options = NosnocProblemOptions();
    solver_options = nosnoc.solver.Options();
    model = NosnocModel();

    problem_options.n_s = 2;
    solver_options.homotopy_update_slope = 0.1;
    problem_options.irk_scheme = irk_scheme;
    
    problem_options.cross_comp_mode = cross_comp_mode;
    problem_options.irk_representation= irk_representation;
    problem_options.dcs_mode = dcs_mode;
    solver_options.print_level = 3;
    solver_options.store_integrator_step_results = 1;
    % discretization parameters
    N_sim = 16;
    T_sim = 1.5;


    problem_options.N_sim = N_sim;
    problem_options.N_finite_elements = 2;
    problem_options.T_sim = T_sim;

    model.x0 = -1;
    x = SX.sym('x',1);
    model.x = x;
    model.c = x;
    model.S = [-1; 1];
    f_1 = [2]; f_2 = [0.2];
    model.F = [f_1 f_2];
    solver_options.use_previous_solution_as_initial_guess = 1;
    integrator = NosnocIntegrator(model, problem_options, solver_options, [], []);
    [results,stats] = integrator.solve();
end
