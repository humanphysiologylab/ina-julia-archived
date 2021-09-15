import BlackBoxOptim: bboptimize, best_candidate

calculate_loss = x -> calculate_loss_segments!(data_buffer, data_true, x, p_kwargs, prob, solve_kwargs_default)
res = bboptimize(calculate_loss; SearchRange=bounds, MaxTime=10. * 60., TraceInterval=10,
                #   PopulationSize=100,
                 );