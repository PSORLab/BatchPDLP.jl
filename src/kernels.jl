
# This kernel generically adds a constraint to each LP, assuming there is one
# line of information in the `constraint_object` for each LP.
function add_LP_constraint_kernel(
    constraint_matrix, 
    right_hand_side, 
    constraint_object, 
    eval_points, 
    active_constraint,
    n_vars::Int32, 
    n_LPs::Int32,
    LP_stride::Int32, 
    index::Int32, 
    geq_flag::Bool, 
    obj_flag::Bool
    )

    # Given a matrix of subgradients of length l, and the points at which the subgradients
    # were calculated, calculate the correct constraints and right-hand sides and add them
    # to the storage at the correct locations (index : LP_stride : end)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Each thread is assigned one constraint line to add (i.e., one LP)
    while idx <= n_LPs
        # If it's the objective function, we need to add a term for the epigraph variable
        # (Note that the objective function should ALWAYS be called as a LEQ constraint)
        if obj_flag
            constraint_matrix[index + LP_stride*(idx - Int32(1)), Int32(1)] = 1.0
        end

        if geq_flag
            # Start with the RHS equal to the negative of the concave relaxation
            RHS_temp = -constraint_object[idx, Int32(2)]

            # For GEQ constraints, we use the concave relaxation subgradient for the
            # constraint matrix, and we add eval_point * CC_subgradient for each dimension
            # to the RHS_temp 
            j = Int32(1)
            while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
                constraint_matrix[index + LP_stride*(idx - Int32(1)), j + Int32(1)] = constraint_object[idx, Int32(4) + (n_vars - Int32(1)) + j]
                RHS_temp += eval_points[idx, j] * constraint_object[idx, Int32(4) + (n_vars - Int32(1)) + j]
                j += Int32(1)
            end
        else
            # Start with the RHS equal to the negative of the convex relaxation 
            # (and then negated again to be positive)
            RHS_temp = constraint_object[idx, Int32(1)]

            # For LEQ constraints, we use the convex relaxation subgradient and negate it for
            # the constraint matrix, and we subtract eval_point * CV_subgradient for each
            # dimension from the RHS_temp 
            j = Int32(1)
            while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
                constraint_matrix[index + LP_stride*(idx - Int32(1)), j + Int32(1)] = -constraint_object[idx, Int32(4) + j]
                RHS_temp -= eval_points[idx, j] * constraint_object[idx, Int32(4) + j]
                j += Int32(1)
            end
        end

        # Finally, set the right-hand side equal to RHS_temp
        right_hand_side[index + LP_stride*(idx - Int32(1))] = RHS_temp

        # Set the active constraint flag to true
        active_constraint[index + LP_stride*(idx - Int32(1))] = true

        idx += stride
    end
    return nothing
end

# This kernel adds a lower bound for the epigraph variable as a constraint,
# assuming there is one line of information in the `obj_result` for each LP.
function add_LP_lower_bound_kernel(
    constraint_matrix, 
    right_hand_side, 
    obj_result, 
    active_constraint, 
    n_LPs::Int32, 
    LP_stride::Int32, 
    index::Int32
    )

    # Given the objective result, which contains [cv, cc, lo, hi, [subgradients...]],
    # pull out the lower bound result and apply it as a constraint on the epigraph variable
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Each thread is assigned one constraint line to add (i.e., one LP)
    while idx <= n_LPs
        # Add the epigraph variable as the only element in that row of the constraint
        constraint_matrix[index + LP_stride*(idx - Int32(1)), Int32(1)] = 1.0

        # The right-hand side value is simply the lower bound from obj_result, since
        # the structure is `constraint_matrix * x >= right_hand_side`
        right_hand_side[index + LP_stride*(idx - Int32(1))] = obj_result[idx, 3]

        # Set the active_constraint row to true
        active_constraint[index + LP_stride*(idx - Int32(1))] = true

        idx += stride
    end
    return nothing
end

# This kernel adds `n_points` constraints to each LP.
# DEPRECATION WARNING: This kernel is outdated and does not change the `active_constraint`
#                      field. This field alerts the main PDLP kernel that the constraint is
#                      in use. If `active_constraint` is not adjusted, the constraint will be
#                      ignored.
# function add_multiple_LP_constraint_kernel(
#     constraint_matrix, 
#     right_hand_side, 
#     constraint_object, 
#     eval_points, 
#     n_vars::Int32, 
#     n_LPs::Int32,
#     LP_stride::Int32, 
#     index::Int32, 
#     geq_flag::Bool, 
#     obj_flag::Bool,
#     n_points::Int32
#     )

#     # Given a matrix of subgradients of length l, and the points at which the subgradients
#     # were calculated, calculate the correct constraints and right-hand sides and add them
#     # to the storage at the correct locations (index : LP_stride : end)
#     idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
#     stride = blockDim().x * gridDim().x

#     # Each thread is assigned one constraint line to add (i.e., one LP)
#     while idx <= n_LPs

#         if geq_flag
#             point = Int32(1)
#             while point <= n_points
#                 constraint_index = (idx - Int32(1))*n_points + point
#                 # Start with the RHS equal to the negative of the concave relaxation
#                 RHS_temp = -constraint_object[constraint_index, Int32(2)]

#                 # For GEQ constraints, we use the concave relaxation subgradient for the
#                 # constraint matrix, and we add eval_point * CC_subgradient for each dimension
#                 # to the RHS_temp 
#                 j = Int32(1)
#                 while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
#                     constraint_matrix[index + LP_stride*(idx - Int32(1)), j + Int32(1)] = constraint_object[constraint_index, Int32(4) + (n_vars - Int32(1)) + j]
#                     RHS_temp += eval_points[constraint_index, j] * constraint_object[constraint_index, Int32(4) + (n_vars - Int32(1)) + j]
#                     j += Int32(1)
#                 end
                
#                 # Finally, set the right-hand side equal to RHS_temp
#                 right_hand_side[index + LP_stride*(idx - Int32(1))] = RHS_temp

#                 index += Int32(1)
#                 point += Int32(1)
#             end
#         else
#             point = Int32(1)
#             while point <= n_points
#                 constraint_index = (idx - Int32(1))*n_points + point
#                 # Start with the RHS equal to the negative of the convex relaxation 
#                 # (and then negated again to be positive)
#                 RHS_temp = constraint_object[constraint_index, Int32(1)]

#                 # If it's the objective function, we need to add a term for the epigraph variable
#                 # (Note that the objective function should ALWAYS be called as a LEQ constraint)
#                 if obj_flag
#                     constraint_matrix[index + LP_stride*(idx - Int32(1)), Int32(1)] = 1.0
#                 end

#                 # For LEQ constraints, we use the convex relaxation subgradient and negate it for
#                 # the constraint matrix, and we subtract eval_point * CV_subgradient for each
#                 # dimension from the RHS_temp 
#                 j = Int32(1)
#                 while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
#                     constraint_matrix[index + LP_stride*(idx - Int32(1)), j + Int32(1)] = -constraint_object[constraint_index, Int32(4) + j]
#                     RHS_temp -= eval_points[constraint_index, j] * constraint_object[constraint_index, Int32(4) + j]
#                     j += Int32(1)
#                 end
                
#                 # Finally, set the right-hand side equal to RHS_temp
#                 right_hand_side[index + LP_stride*(idx - Int32(1))] = RHS_temp
            
#                 index += Int32(1)
#                 point += Int32(1)
#             end
#         end

#         idx += stride
#     end
#     return nothing
# end

# This kernel recognizes that `n_points` relaxations exist for each LP, and
# then selects only the first point to source lower bound information, since
# the information should be the same for every evaluation within one LP.
function add_multiple_LP_lower_bound_kernel(
    constraint_matrix,
    right_hand_side,
    obj_result,
    n_LPs::Int32,
    LP_stride::Int32,
    index::Int32,
    active_constraint,
    n_points::Int32
    )

    # Given the objective result, which contains [cv, cc, lo, hi, [subgradients...]],
    # pull out the lower bound result and apply it as a constraint on the epigraph variable
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Each thread is assigned one constraint line to add (i.e., one LP)
    while idx <= n_LPs
        # Add the epigraph variable as the only element in that row of the constraint
        constraint_matrix[index + LP_stride*(idx - Int32(1)), Int32(1)] = 1.0

        # The right-hand side value is simply the lower bound from obj_result, since
        # the structure is `constraint_matrix * x >= right_hand_side`. Since we have
        # n_points per LP, take only the bound from the first of the points
        right_hand_side[index + LP_stride*(idx - Int32(1))] = obj_result[Int32((idx - Int32(1))*n_points) + Int32(1), 3]

        # Set the active_constraint row to true
        active_constraint[index + LP_stride*(idx - Int32(1))] = true

        idx += stride
    end
    return nothing
end

# This kernel scans through `n_points` relaxations for each LP to identify the "best"
# constraint, before adding it as a constraint in the LP. It does this process
# `num_linearizations` times. Note that if no more unique constraints can be added
# (or if existing constraints would only provide similar information), the constraint
# matrix will be filled in with 0s and `active_constraint` will be marked as `false`.
function add_best_obj_LP_constraints_kernel(
    constraint_matrix, 
    right_hand_side, 
    constraint_object, 
    eval_points, 
    comparison_vector,
    subgradient_checksum,
    n_vars::Int32, 
    n_LPs::Int32,
    LP_stride::Int32, # Total length of each LP
    LP_index::Int32, # Index of the most recent constraint
    active_constraint,
    n_points::Int32, # Total number of linearizations to look through
    num_linearizations::Int32 # Number of linearizations to use
    )

    # Given a matrix of subgradients of length l, and the points at which the subgradients
    # were calculated, calculate the correct constraints and right-hand sides and add them
    # to the storage at the correct locations (index : LP_stride : end)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Each thread is assigned to one LP
    while idx <= n_LPs
        # Loop through num_linearizations times before moving on
        lin_ID = Int32(1)
        while lin_ID <= num_linearizations
            # First, find the lowest value and its index in `comparison_vector`.
            # Note that the vector is a sum of squared values, so it's always positive.
            # I.e., this picks the subgradient closest to 0.0 (not the most negative).
            min_value = Inf
            min_index = Int32(0)
            point = Int32(1)
            while point <= n_points
                val = comparison_vector[(idx - Int32(1))*n_points + point]
                if val < min_value
                    # If the value is an improvement, set it as the new best
                    min_value = val
                    min_index = (idx - Int32(1))*n_points + point
                end
                # (If it's a tie, the choice is arbitrary, so we won't address that case)
                point += Int32(1)
            end

            # If we got a minimum value that's not Inf, then we should add an LP constraint.
            if ~isinf(min_value)
                # Start with the RHS equal to the negative of the convex relaxation 
                # (and then negated again to be positive)
                RHS_temp = constraint_object[min_index, Int32(1)]

                # Since this is for the objective function, we need to add a term for the epigraph variable
                constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), Int32(1)] = 1.0

                # For LEQ constraints, we use the convex relaxation subgradient and negate it for
                # the constraint matrix, and we subtract eval_point * CV_subgradient for each
                # dimension from the RHS_temp 
                j = Int32(1)
                while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
                    constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), j + Int32(1)] = -constraint_object[min_index, Int32(4) + j]
                    RHS_temp -= eval_points[min_index, j] * constraint_object[min_index, Int32(4) + j]
                    j += Int32(1)
                end
                
                # Finally, set the right-hand side equal to RHS_temp
                right_hand_side[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = RHS_temp

                # Set the active constraint flag to true
                active_constraint[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = true

            # If we got a minimum value of Inf, then we have no more LP constraints to add. We'll 
            # put in all zeros instead.
            else
                # We can immediately go to the right-hand side and constraint matrix and set the values to 0
                j = Int32(1)
                constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), Int32(1)] = 0.0
                while j < n_vars
                    constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), j + Int32(1)] = 0.0
                    j += Int32(1)
                end
                right_hand_side[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = 0.0

                # Set the active constraint flag to false
                active_constraint[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = false

                # Note that we can't skip the rest of the lin_ID's here, since future lines in the constraint matrix
                # may need to be reset to 0 from previous calculations.
            end

            # Now we must clean up the comparison vector so that the same subgradient doesn't get used
            # again. We will do this by checking the subgradient_checksum. If a value in that vector
            # matches the value at min_index, we set comparison_vector at that index to Inf.
            point = Int32(1)
            min_checksum = subgradient_checksum[min_index]
            while point <= n_points
                if subgradient_checksum[(idx - Int32(1))*n_points + point] == min_checksum
                    comparison_vector[(idx - Int32(1))*n_points + point] = Inf
                end
                point += Int32(1)
            end

            # ALTERNATIVE FOR TESTING PURPOSES ONLY: 
            # comparison_vector[min_index] = Inf

            lin_ID += Int32(1)
        end
        idx += stride
    end
    return nothing
end

# This kernel scans through `n_points` relaxations for each LP to identify the "best"
# constraint, before adding it as a constraint in the LP. It does this process
# `num_linearizations` times. Note that if no more unique constraints can be added
# (or if existing constraints would only provide similar information), the constraint
# matrix will be filled in with 0s and `active_constraint` will be marked as `false`.
function add_best_cons_LP_constraints_kernel(
    constraint_matrix, 
    right_hand_side, 
    constraint_object, 
    eval_points, 
    comparison_vector,
    subgradient_checksum,
    n_vars::Int32, 
    n_LPs::Int32,
    LP_stride::Int32, # Total length of each LP
    LP_index::Int32, # Index of the most recent constraint
    active_constraint,
    geq_flag::Bool, # Indicator for "GEQ" constraint in original NLP
    n_points::Int32, # Total number of linearizations to look through
    num_linearizations::Int32 # Number of linearizations to use
    )

    # Given a matrix of subgradients of length l, and the points at which the subgradients
    # were calculated, calculate the correct constraints and right-hand sides and add them
    # to the storage at the correct locations (index : LP_stride : end)
    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Each thread is assigned to one LP
    while idx <= n_LPs
        # Loop through num_linearizations times before moving on
        lin_ID = Int32(1)
        while lin_ID <= num_linearizations
            # First, find the lowest absolute value and its index in `comparison_vector`.
            # This vector contains either the convex or concave relaxation value, and we
            # want the value closest to 0, hence why we use the absolute value.
            min_value = Inf
            min_index = Int32(0)
            point = Int32(1)
            while point <= n_points
                val = abs(comparison_vector[(idx - Int32(1))*n_points + point])
                if val < min_value
                    # If the value is an improvement, set it as the new best
                    min_value = val
                    min_index = (idx - Int32(1))*n_points + point
                end
                point += Int32(1)
            end

            # If we got a value that isn't Inf, then we have a constraint to add
            if ~isinf(min_value)
                # Apply the constraint associated with `min_index`
                if geq_flag
                    # Start with the RHS equal to the negative of the concave relaxation
                    RHS_temp = -constraint_object[min_index, Int32(2)]

                    # Since this is for a constraint, make sure the entry in the
                    # constraint matrix for the epigraph variable is 0.
                    constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), Int32(1)] = 0.0

                    # For GEQ constraints, we use the concave relaxation subgradient for the
                    # constraint matrix, and we add eval_point * CC_subgradient for each dimension
                    # to the RHS_temp
                    j = Int32(1)
                    while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
                        constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), j + Int32(1)] = constraint_object[min_index, Int32(4) + (n_vars - Int32(1)) + j]
                        RHS_temp += eval_points[min_index, j] * constraint_object[min_index, Int32(4) + (n_vars - Int32(1)) + j]
                        j += Int32(1)
                    end
                    
                    # Finally, set the right-hand side equal to RHS_temp
                    right_hand_side[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = RHS_temp

                    # Set the active constraint flag to true
                    active_constraint[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = true
                else
                    # Start with the RHS equal to the negative of the convex relaxation 
                    # (and then negated again to be positive)
                    RHS_temp = constraint_object[min_index, Int32(1)]

                    # Since this is for a constraint, make sure the entry in the
                    # constraint matrix for the epigraph variable is 0.
                    constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), Int32(1)] = 0.0

                    # For LEQ constraints, we use the convex relaxation subgradient and negate it for
                    # the constraint matrix, and we subtract eval_point * CV_subgradient for each
                    # dimension from the RHS_temp 
                    j = Int32(1)
                    while j < n_vars # Note: <, not <=, because n_vars includes the epigraph variable
                        constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), j + Int32(1)] = -constraint_object[min_index, Int32(4) + j]
                        RHS_temp -= eval_points[min_index, j] * constraint_object[min_index, Int32(4) + j]
                        j += Int32(1)
                    end
                    
                    # Finally, set the right-hand side equal to RHS_temp
                    right_hand_side[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = RHS_temp

                    # Set the active constraint flag to true
                    active_constraint[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = true
                end

            # If we got a minimum value of Inf, then we have no more LP constraints to add. We'll 
            # put in all zeros instead.
            else
                # We can immediately go to the right-hand side and constraint matrix and set the values to 0
                j = Int32(1)
                constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), Int32(1)] = 0.0
                while j < n_vars
                    constraint_matrix[LP_index + lin_ID + LP_stride*(idx - Int32(1)), j + Int32(1)] = 0.0
                    j += Int32(1)
                end
                right_hand_side[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = 0.0

                # Set the active constraint flag to false
                active_constraint[LP_index + lin_ID + LP_stride*(idx - Int32(1))] = false

                # Note that we can't skip the rest of the lin_ID's here, since future lines in the constraint matrix
                # may need to be reset to 0 from previous calculations.
            end

            # Now we must clean up the comparison vector so that the same subgradient doesn't get used
            # again. We will do this by checking the subgradient_checksum. If a value in that vector
            # matches the value at min_index, we set comparison_vector at that index to Inf.
            point = Int32(1)
            min_checksum = subgradient_checksum[min_index]
            while point <= n_points
                if subgradient_checksum[(idx - Int32(1))*n_points + point] == min_checksum
                    comparison_vector[(idx - Int32(1))*n_points + point] = Inf
                end
                point += Int32(1)
            end
            
            lin_ID += Int32(1)
        end
        idx += stride
    end
    return nothing
end


# Several kernels for calculating overall rescaling prior to running PDLP. Calculations
# are performed separately for each LP.
function ruiz_variable_kernel(
    result_storage, 
    constraint_matrix, 
    current_LP_length, 
    total_LP_length
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Get size information from the constraint matrix
    len = Int32(size(constraint_matrix, 1))
    width = Int32(size(constraint_matrix, 2))

    # Use the constraint matrix size, and the size of each LP, to figure out
    # how high the indexing should go
    maxlen = Int32((len*width)/total_LP_length)

    # For each column of each LP, get the maximum value of the absolute values of
    # the elements in that column, and set the result storage for that LP and
    # variable to the sqrt of that value
    @inbounds while idx <= maxlen
        # Figure out which LP is being worked on
        LP = Int32(cld(idx,width))

        # Initialize the row/col information, with row starting at 2 since we fill
        # in `maxval` before the loop.
        row = Int32(2)
        col = idx - Int32((LP-Int32(1))*width)

        # Calculate the maximum value in each column, for each LP
        maxval = abs(constraint_matrix[(LP - Int32(1))*total_LP_length + Int32(1), col])
        while row <= current_LP_length
            maxval = max(abs(constraint_matrix[(LP - Int32(1))*total_LP_length + row, col]), maxval)
            row += Int32(1)
        end

        # Set the value to 1.0 if the maximum was zero.
        if maxval < 1e-12
            maxval = 1.0
        end

        # Store the results
        result_storage[LP,col] = sqrt(maxval)
        idx += stride
    end
    return nothing
end

function ruiz_constraint_kernel(
    result_storage, 
    constraint_matrix, 
    current_LP_length, 
    total_LP_length
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x
    len = Int32(size(constraint_matrix, 1)) # number of constraints
    width = Int32(size(constraint_matrix, 2)) # number of variables

    # For each row of the constraint matrix, get the maximum value of the absolute 
    # values of the elements in that row, and set the result storage to the sqrt of 
    # that value
    @inbounds while idx <= len
        # Skip rows that aren't currently active.
        new_ID = idx + (cld(idx,current_LP_length)-1)*(total_LP_length - current_LP_length)
        if new_ID > len
            break
        end

        col = Int32(1)
        maxval = 0.0
        while col <= width
            maxval = max(abs(constraint_matrix[new_ID, col]), maxval)
            col += Int32(1)
        end
        if maxval < 1e-12
            maxval = 1.0
        end
        result_storage[new_ID] = sqrt(maxval)
        idx += stride
    end
    return nothing
end

function scaling_part1_kernel(
    objective_vector, 
    variable_lower_bounds, 
    variable_upper_bounds, 
    variable_rescaling
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    len = size(objective_vector, 1)
    width = size(objective_vector, 2)

    @inbounds while idx <= len
        var = Int32(1)
        while var <= width
            objective_vector[idx, var] /= variable_rescaling[idx, var]
            variable_lower_bounds[idx, var] *= variable_rescaling[idx, var]
            variable_upper_bounds[idx, var] *= variable_rescaling[idx, var]
            var += Int32(1)
        end
        idx += stride
    end
    return nothing
end

function scaling_part2_kernel(
    right_hand_side, 
    constraint_matrix, 
    current_LP_length, 
    total_LP_length, 
    variable_rescaling, 
    constraint_rescaling
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    len = size(constraint_matrix, 1)
    width = size(constraint_matrix, 2)

    @inbounds while idx <= len
        # Skip rows that aren't currently active.
        new_ID = idx + (cld(idx,current_LP_length)-1)*(total_LP_length - current_LP_length)
        if new_ID > len
            break
        end

        right_hand_side[new_ID] /= constraint_rescaling[new_ID]
        LP = Int32(cld(new_ID,total_LP_length))
        var = Int32(1)
        while var <= width
            # Each ROW of the constraint matrix is divided by that row's constraint rescaling
            constraint_matrix[new_ID,var] /= constraint_rescaling[new_ID]

            # Each COL of the constraint matrix is divided by that variable's variable rescaling
            constraint_matrix[new_ID,var] /= variable_rescaling[LP,var]
            var += Int32(1)
        end
        idx += stride
    end
    return nothing
end

function pock_chambolle_variable_kernel(
    result_storage, 
    constraint_matrix, 
    current_LP_length, 
    total_LP_length, 
    alpha
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    # Get size information from the constraint matrix
    len = Int32(size(constraint_matrix, 1))
    width = Int32(size(constraint_matrix, 2))

    # Use the constraint matrix size, and the size of each LP, to figure out
    # how high the indexing should go
    maxlen = Int32((len*width)/total_LP_length)

    # For each column of the LP, sum up abs(x)^(2 - alpha), and save sqrt(result)
    # to result_storage
    @inbounds while idx <= maxlen
        # Figure out which LP is being worked on
        LP = Int32(cld(idx,width))

        # Initialize the row/col information, with row starting at 2 since we fill
        # in `result` before the loop.
        row = Int32(2)
        col = idx - Int32((LP-Int32(1))*width)

        # Calculate the result in each column, for each LP
        result = abs(constraint_matrix[(LP - Int32(1))*total_LP_length + Int32(1), col])^(2 - alpha)
        while row <= current_LP_length
            result += (abs(constraint_matrix[(LP - Int32(1))*total_LP_length + row, col]))^(2 - alpha)
            row += Int32(1)
        end

        # Set the value to 1.0 if the result was zero.
        if result < 1e-12
            result = 1.0
        end

        # Store the results
        result_storage[LP,col] = sqrt(result)
        idx += stride
    end
    return nothing
end

function pock_chambolle_constraint_kernel(
    result_storage, 
    constraint_matrix, 
    current_LP_length, 
    total_LP_length,
    alpha, 
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x
    len = Int32(size(constraint_matrix, 1))
    width = Int32(size(constraint_matrix, 2))

    # For each row of the constraint matrix, sum up abs(x)^alpha, and save
    # sqrt(result) to result_storage
    @inbounds while idx <= len
        # Skip rows that aren't currently active.
        new_ID = idx + (cld(idx,current_LP_length)-1)*(total_LP_length - current_LP_length)
        if new_ID > len
            break
        end

        col = Int32(2)
        result = abs(constraint_matrix[new_ID, Int32(1)])^alpha
        while col <= width
            result += abs(constraint_matrix[new_ID, col])^alpha
            col += Int32(1)
        end
        if result < 1e-12
            result = 1.0
        end
        result_storage[new_ID] = sqrt(result)
        idx += stride
    end
    return nothing
end

# This kernel calculates the initial primal weight by taking the norm of the objective
# vector and the right-hand side, for each LP, and using their ratio to set the
# primal weight as cupdlp.jl does
function primal_weight_kernel(
    result, 
    objective_vector, 
    right_hand_side, 
    n_groups, 
    group_length, 
    current_length
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x
    width = Int32(size(objective_vector, 2))

    while idx <= n_groups
        # Calculate the objective norm first
        obj_norm = 0.0
        col_or_row = Int32(1)
        while col_or_row <= width
            obj_norm += objective_vector[idx,col_or_row]^2
            col_or_row += Int32(1)
        end

        # Then calculate the right-hand side norm
        rhs_norm = 0.0
        col_or_row = Int32(1)
        while col_or_row <= current_length
            rhs_norm += right_hand_side[(idx - Int32(1))*group_length + col_or_row]^2
            col_or_row += Int32(1)
        end

        # And finally, calculate the primal importance and save it to the result
        if obj_norm > 0.0 && rhs_norm > 0.0
            result[idx] = sqrt(obj_norm)/sqrt(rhs_norm)
        else
            result[idx] = 1.0
        end
        idx += stride
    end
    return nothing
end

# This kernel calculates the initial primal weight by taking the norm of the objective
# vector and the right-hand side, for each LP, and using their ratio to set the
# primal weight as cupdlpx does
function primal_weight_kernel_rHalpern(
    result, 
    objective_vector, 
    right_hand_side, 
    n_groups, 
    group_length, 
    current_length
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x
    width = Int32(size(objective_vector, 2))

    while idx <= n_groups
        # Calculate the objective norm first
        obj_norm = 0.0
        col_or_row = Int32(1)
        while col_or_row <= width
            obj_norm += objective_vector[idx,col_or_row]^2
            col_or_row += Int32(1)
        end

        # Then calculate the right-hand side norm
        rhs_norm = 0.0
        col_or_row = Int32(1)
        while col_or_row <= current_length
            rhs_norm += right_hand_side[(idx - Int32(1))*group_length + col_or_row]^2
            col_or_row += Int32(1)
        end

        # And finally, calculate the primal importance and save it to the result
        if obj_norm > 0.0 && rhs_norm > 0.0
            # result[idx] = sqrt(obj_norm)/sqrt(rhs_norm)
            result[idx] = (sqrt(obj_norm) + 1) / (sqrt(rhs_norm) + 1)
        else
            result[idx] = 1.0
        end
        
        idx += stride
    end
    return nothing
end

# This function finds 1/norm(matrix, Inf) by scanning each row up to current_length
# in each LP for the maximum.
function group_max_kernel(
    result, 
    matrix, 
    n_LPs, 
    total_LP_length, 
    current_LP_length
    )

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x
    len = size(matrix, 1)
    width = size(matrix, 2)

    while idx <= n_LPs
        row = Int32(1)
        maximum = -Inf
        while row <= current_LP_length
            col = Int32(1)
            while col <= width
                maximum = max(abs(matrix[(idx - Int32(1))*total_LP_length + row, col]), maximum)
                col += Int32(1)
            end
            row += Int32(1)
        end
        result[idx] = 1.0/maximum
        idx += stride
    end
    return nothing
end

# this kernel function is for doing power iterations to estimate ||G||_2. 
function group_power_kernel(
    result, # [n_LPs] - output step sizes
    matrix, 
    n_LPs, 
    n_vars,
    nz_count,
    nz_rows,
    nz_cols,
    active_constraint,
    total_LP_length, 
    current_LP_length,
    u_vec, # [n_LPs, n_vars] - intermediate: A·v
    eigenvector_d, # [n_LPs * total_LP_length] - current eigenvector
    new_eigenvector_d, # [n_LPs * total_LP_length] - next eigenvector (A·A^T·v)
    iterations,
    tolerance
    )

    LP = blockIdx().x
    idx = threadIdx().x

    block_stride = blockDim().x
    grid_stride = gridDim().x

    shared_space = @cuDynamicSharedMem(Float64, max(current_LP_length, n_vars))


    var_stride = Int32(1) << floor(Int32, log2(n_vars))
    len_stride = Int32(1) << floor(Int32, log2(current_LP_length))

    while LP <= n_LPs
        
        # LP-specific starting position in the 1D workspace array (1-based)
        active_row = (LP - Int32(1)) * total_LP_length

        iter = Int32(1)
        sigma_max_sq = 1.0
        
        # Power iterations
        while iter <= iterations 

            #########################################################
            # STEP 1: Copy and normalize guess_vector to new_vector #
            #########################################################
            # Copy eigenvector_d into new_eigenvector_d
            while idx <= current_LP_length
                new_eigenvector_d[active_row + idx] = eigenvector_d[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            sync_threads()

            # Compute L2 norm of new_eigenvector_d: ||new_eigenvector_d||_2

            while idx <= current_LP_length
                shared_space[idx] = new_eigenvector_d[active_row + idx] ^ 2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            
            eigenvector_norm = sqrt(shared_space[1])
            inv_eigenvector_norm = 1 / eigenvector_norm
            # Normalize: new_vector = new_vector / ||new_vector||_2

            while idx <= current_LP_length
                new_eigenvector_d[active_row + idx] =  new_eigenvector_d[active_row + idx] * inv_eigenvector_norm
                idx += block_stride
            end

            idx = threadIdx().x
            sync_threads()

            ##############################################################################
            # STEP 2: Compute u_vec = A^T * new_eigenvector_d ([n*1] = [n*m] x [m*1])    #
            ##############################################################################

            while idx <= n_vars 
                shared_space[idx] = 0.0
                idx += block_stride
            end
            idx = threadIdx().x
            while idx <= nz_count
                if active_constraint[active_row + nz_rows[idx]]
                    CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), matrix[active_row + nz_rows[idx], nz_cols[idx]] * new_eigenvector_d[active_row + nz_rows[idx]])
                end
                idx += block_stride
            end
            idx = threadIdx().x
            sync_threads()

            while idx <= n_vars
                u_vec[LP, idx] = shared_space[idx]
                idx += block_stride
            end
            idx = threadIdx().x
            sync_threads()

            ################################################################################
            # STEP 3: Compute eigenvector_d = A * u_vec ([m*1] = [m*n] x [n*1])            #
            ################################################################################

            while idx <= current_LP_length 
                shared_space[idx] = 0.0
                idx += block_stride
            end
            idx = threadIdx().x
            while idx <= nz_count
                if active_constraint[active_row + nz_rows[idx]]
                    CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), matrix[active_row + nz_rows[idx], nz_cols[idx]] * u_vec[LP, nz_cols[idx]])
                end
                idx += block_stride
            end
            idx = threadIdx().x
            sync_threads()

            while idx <= current_LP_length
                eigenvector_d[active_row + idx] = shared_space[idx]
                idx += block_stride
            end
            
            idx = threadIdx().x
            sync_threads()

            #########################################################
            # STEP 4: Compute Rayleigh quotient                     #
            # sigma_max_sq = guess_vector * next_vector             #
            #########################################################

            while idx <= current_LP_length
                shared_space[idx] = new_eigenvector_d[active_row + idx] * eigenvector_d[active_row + idx]
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            
            sigma_max_sq = shared_space[1]
            sync_threads()

            #########################################################
            # STEP 5: Compute residual for convergence check        #
            # residual = guess_vector - sigma_max_sq * new_vector   #
            #########################################################

            while idx <= current_LP_length
                shared_space[idx] = (eigenvector_d[active_row + idx] - sigma_max_sq * new_eigenvector_d[active_row + idx]) ^ 2
                idx += block_stride
            end
            idx = threadIdx().x
            parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
            sync_threads()

            residual_norm = sqrt(shared_space[1])

            # Check convergence
            if residual_norm < tolerance
                break
            end

            iter += Int32(1)
        end
        
        # Compute step size: tau = 0.998 / sqrt(sigma_max_sq)
        if idx == 1
            result[LP] = 0.998 / sqrt(sigma_max_sq)
            
        end
        LP += grid_stride
    end
    return nothing
end

function compute_bound_contrib_kernel(
    results,
    right_hand_side,
    total_LP_length,
)

    idx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    stride = blockDim().x * gridDim().x

    while idx <= total_LP_length
        Li = right_hand_side[idx]
        acc = 0.0 
        
        # Since there are no upper bounds, we only check if the lower bound (RHS) is finite.
        if isfinite(Li)
            acc += Li * Li
        end
        
        results[idx] = acc
        
        # CRITICAL: Advance the index to prevent an infinite loop
        idx += stride
    end
    
    return nothing 
end