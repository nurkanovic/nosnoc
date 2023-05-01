function print_problem_details(results,model,solver_initialization,filename)
if exist('filename') && ~isempty(filename)
    delete(filename);
    fileID = fopen(filename, 'w');
else
    fileID = 1;
end
fprintf(fileID, "i\t\t lbg\t\t\t ubg \t\t\t g_val \t\t infeasible \t\t g_exp\n");
inf_trashhold = 1e-3;
for i = 1:length(solver_initialization.lbg)
    expr_str = formattedDisplayText(model.g(i));
    expr_val = full(results.g(i));
    if expr_val < solver_initialization.lbg(i)-inf_trashhold  || expr_val > solver_initialization.ubg(i)+inf_trashhold
        str_inf = '[x]';
    else
        str_inf = '  ';
    end
    fprintf(fileID, "%d\t\t %6.2e \t\t %6.2e \t\t %6.2e \t\t  %s \t\t  %s \n", i, solver_initialization.lbg(i), solver_initialization.ubg(i), expr_val,str_inf, expr_str);
end

fprintf(fileID, "i\t\t lbw\t\t\t ubw \t\t\t w_val \t\t infeasible \t\t w\n");
inf_trashhold = 1e-3;
for i = 1:length(solver_initialization.lbw)
    expr_str = formattedDisplayText(model.w(i));
    expr_val = full(results.x(i));
    if expr_val < solver_initialization.lbw(i)-inf_trashhold  || expr_val > solver_initialization.ubw(i)+inf_trashhold
        str_inf = '[x]';
    else
        str_inf = '  ';
    end
    fprintf(fileID, "%d\t\t %6.2e \t\t %6.2e \t\t %6.2e \t\t  %s \t\t  %s \n", i, solver_initialization.lbw(i), solver_initialization.ubw(i), expr_val,str_inf, expr_str);
end

% fprintf(fileID, "\nw\t\t\tw0\t\tlbw\t\tubw\n");
% for i = 1:length(obj.lbw)
%     % keyboard
%     expr_str = pad(formattedDisplayText(obj.w(i)), 20);
%     lb_str = pad(sprintf('%.2e', obj.lbw(i)), 10);
%     fprintf(fileID, "%s\t%.2e\t%s\t%.2e\t\n", expr_str, obj.w0(i), lb_str, obj.ubw(i));
% end

% fprintf(fileID, '\nobjective\t\t objective_val\n');
% fprintf(fileID, strcat(formattedDisplayText(obj.cost), '\n'));
end