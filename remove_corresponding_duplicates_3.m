 %Specify the filename and sheet name
input_filename = 'side_collision_same_direction_all_files_approaching angles';
filename = append(input_filename,'.xlsx');
toatal_sheets =27;
sheet_no = 1;

while sheet_no <= toatal_sheets
    % Read the data from the first four columns (A, B, C, and D)
    [num,txt,raw] = xlsread(filename, sheet_no, 'A:D');

    % Assign the data to variables with the specified names
    x_traj1 = num(:,1);
    y_traj1 = num(:,2);
    x_traj2 = num(:,3);
    y_traj2 = num(:,4);

    % Round off input arrays
    x_traj1 = round(x_traj1, 2);
    y_traj1 = round(y_traj1, 2);
    x_traj2 = round(x_traj2, 2);
    y_traj2 = round(y_traj2, 2);

    % Remove duplicate coordinates in traj1 and corresponding coordinates in traj2
    [x_traj1_intermediate, idx1] = unique(x_traj1, 'stable');
    y_traj1_intermediate = y_traj1(idx1);
    x_traj2_intermediate = x_traj2(idx1);
    y_traj2_intermediate = y_traj2(idx1);

    % Remove duplicate coordinates in traj2_intermediate and corresponding coordinates in traj1_intermediate
    [x_traj2_final, idx2] = unique(x_traj2_intermediate, 'stable');
    x_traj1_final = x_traj1_intermediate(idx2);
    y_traj1_final = y_traj1_intermediate(idx2);
    y_traj2_final = y_traj2_intermediate(idx2);

    % Display final trajectories
    % disp("Final Trajectory 1:");
    % disp([x_traj1_final; y_traj1_final]);
    % disp("Final Trajectory 2:");
    % disp([x_traj2_final; y_traj2_final])


    % Determine whether to take increasing or decreasing values of x
    if x_traj1_final(end) - x_traj1_final(1) > 0
        inc_order = true;
    else
        inc_order = false;
    end

    x_traj1_final_inter = x_traj1_final;
    x_traj2_final_inter = x_traj2_final;
    y_traj1_final_inter = y_traj1_final;
    y_traj2_final_inter = y_traj2_final;

    while 1
        % Remove indices from trajectory 1
        if inc_order
            indices_to_remove = find(diff(x_traj1_final_inter) < 0) + 1;
        else
            indices_to_remove = find(diff(x_traj1_final_inter) > 0) + 1;
        end

%         disp(length(indices_to_remove));
        if isempty(indices_to_remove) 
            break;
        end

        x_traj1_final_inter(indices_to_remove) = [];
        y_traj1_final_inter(indices_to_remove) = [];
        
        % Remove corresponding indices from trajectory 2
        x_traj2_final_inter(indices_to_remove) = [];
        y_traj2_final_inter(indices_to_remove) = [];
    end

%     disp('loop-2');
    % Determine whether to take increasing or decreasing values of x
    if x_traj2_final_inter(end) - x_traj2_final_inter(1) > 0
        inc_order = true;
    else
        inc_order = false;
    end

    x_traj2_final_final = x_traj2_final_inter;
    y_traj2_final_final = y_traj2_final_inter;
    x_traj1_final_final = x_traj1_final_inter;
    y_traj1_final_final = y_traj1_final_inter;

    while 1
        % Remove indices from trajectory 2
        if inc_order
            indices_to_remove = find(diff(x_traj2_final_final) < 0) + 1;
        else
            indices_to_remove = find(diff(x_traj2_final_final) > 0) + 1;
        end
%         disp(length(indices_to_remove));
        if isempty(indices_to_remove) 
            break;
        end


        x_traj2_final_final(indices_to_remove) = [];
        y_traj2_final_final(indices_to_remove) = [];

        % Remove corresponding indices from trajectory 1
        x_traj1_final_final(indices_to_remove) = [];
        y_traj1_final_final(indices_to_remove) = [];
    end
%     disp('loop-2 ended');
    
    % fit a quadrtic equation in both trajectories
    [y1_hat,gof1] = fit(x_traj1_final_final,y_traj1_final_final,'poly2');
    [y2_hat,gof2] = fit(x_traj2_final_final,y_traj2_final_final,'poly2');
    
%     disp(y1_hat);
%     disp(class(y1_hat));
%     disp(class(x_traj1_final_final));
    % create a matrix with four columns
    
    y1_hat_final = feval(y1_hat, x_traj1_final_final);
    y2_hat_final = feval(y2_hat, x_traj2_final_final);

    data = [x_traj1_final_final, y_traj1_final_final, x_traj2_final_final, y_traj2_final_final,y1_hat_final,y2_hat_final];
    % get the size of the matrix
    [rows, cols] = size(data);

    % specify the range of cells
    range = sprintf('A1:F%d', rows);

    % write the matrix to an Excel file
    output_filename = append(input_filename,'_output');
    output_sheetname = append('sheet',string(sheet_no));
    xlswrite(output_filename, data, output_sheetname, range);
    
    % printing details of curve fitting 
    fprintf('sheet %d completed.\n',sheet_no);
    disp('trajectory 1');
    disp(gof1);
    disp('trajectory 2');
    disp(gof2);
sheet_no = sheet_no +1;
end
