import json
import copy
import logging
import argparse
import os
import sys

from time import time
from shutil import copy
import gurobipy as grb

from datetime import datetime

now = datetime.now()
try:
    from MIPLIBing import MIPLIBing
    from MIPLIBing import Libraries
    from MIPLIBing import Status
except:
    pass
import global_support

def dump_if_over_time(master, certificate,
                      reduced_certificate, reduced,
                      target_directory,name_of_instance):
    if time()- global_support.STARTING_TIME>= global_support.TIME_LIMIT:
        logging.debug("too slow")
        certificate["finished"] = "time_kill"
        with open(os.path.join(target_directory,
                               name_of_instance.split(".")[0] + ".certificate"), 'w') as output_file:
            output_file.write(json.dumps(certificate))
        exit()


def downloader(name_of_instance):
    """Download instance automatically from web. MIPLib 2017."""
    logging.debug(f'Downloading instance {name_of_instance}')
    logging.debug(name_of_instance)
    mip = MIPLIBing(library = Libraries.MIPLIB2017_Collection)
    # mip = MIPLIBing(Libraries.MINLPLIB, file_extension = "lp")
    instance = mip.get_instances(name_of_instance)


def cb_tracking(model, where, time_list):
    """Less important time tracking callback."""
    if where == grb.GRB.Callback.MESSAGE:
        if not hasattr(model, '_last_memory_check'):
            model._last_memory_check=time()
        if model._last_memory_check - time()<=-5:
            logging.debug(f"Memory_usage_currently_{float(os.popen('free -t -m').readlines()[-1].split()[1:][1])}")
            if float(os.popen('free -t -m').readlines()[-1].split()[1:][1]) >= global_support.MEMORY_LIMIT:
                logging.debug("too large")
                model.terminate()
            model._last_memory_check=time()
    if where == grb.GRB.Callback.MIPSOL:
        time_list.append(time())


def cb_subsolving(model, where,
                  feasible_sub,
                  infeasible_sub,
                  certificate,
                  expr_terms,target_directory, name_of_instance):
    """
    Sub Problem Callback.

    Add one certificate constraint to master model whenever a feasible
    solution with positive value is found.

    Parameters
    ----------
    model : gurobipy.Model
        master model.
    where : gurobipy Callback Code
        see gurobi doc
    feasible_sub : gurobipy.Model
        the feasible-sub problem
    infeasible_sub : gurobipy.Model
        the infeasible-sub problem
    certificate : dict
        collects certificate information
    expr_terms : dict
        Condensed information about certificate interpreted as constraints.

    Returns
    -------
    None.

    """
    dump_if_over_time(model, certificate,
                      "Foo", False,
                      target_directory, name_of_instance)
    if where == grb.GRB.Callback.MESSAGE:
        if not hasattr(model, '_last_memory_check'):
            model._last_memory_check=time()
        if model._last_memory_check - time()<=-30:
            logging.debug(f"Memory_usage_currently_{float(os.popen('free -t -m').readlines()[-1].split()[1:][1])}")
            if float(os.popen('free -t -m').readlines()[-1].split()[1:][1]) >= global_support.MEMORY_LIMIT:
                logging.debug("too large")
                certificate["finished"] = "memory_kill"
                model.terminate()
            model._last_memory_check=time()
    if where == grb.GRB.Callback.MIPSOL:
        logging.debug(f"Time_usage_currently_{time()-global_support.STARTING_TIME}")
        dump_if_over_time(model, certificate,
                          "Foo", False,target_directory, name_of_instance)

        logging.debug('Entering callback function')
        if model.cbGetSolution(model._epsilon) <= -global_support.SAFETY_PARAM*int(cb_subsolving.tmp_calls>5)*max(model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND), - model._epsilon.ub):
            logging.info("Iteration: " + str(cb_subsolving.calls) + ", Current Progress: "
                         + str(model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND))
                         + ", Epsilon: " + str(model.cbGetSolution(model._epsilon))
                         + "invalidinvalidinvalidinvalidinvalid")
            logging.debug('Epsilon is equal/close to 0. Leaving Callback.')
            if model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND) >= - 1e-06:
                certificate["finished"] = "regular"
            return
        if cb_subsolving.tmp_calls == global_support.RESET:
            cb_subsolving.tmp_calls = 0
            model.terminate()
        logging.info("Iteration: " + str(cb_subsolving.calls) + ", Current Progress: "
                     + str(model.cbGet(grb.GRB.Callback.MIPSOL_OBJBND))
                     + ", Epsilon: " + str(model.cbGetSolution(model._epsilon)))
        logging.debug(f'Master objective value: {- model.cbGet(grb.GRB.Callback.MIPSOL_OBJ)}')
        iteration = cb_subsolving.calls
        tmp_cert = certificate['iterations']['iteration' + str(iteration)] = {"nbr": iteration}
        tmp_cert["begin_time"] = time() - global_support.STARTING_TIME
        tmp_cert['suspension point'] = {variable.VarName: model.cbGetSolution(variable)
                                            for variable in model._vars
                                            if model.cbGetSolution(variable) != 0}
        for i in tmp_cert.get('suspension point'):
            logging.debug(i)
        feasible_sub.update()
        for variable in model._vars:
            feasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                round(model.cbGetSolution(variable))
        feasible_sub.update()
        tmp_cert["feas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"]
        feasible_sub.optimize()
        tmp_cert["feas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["feas_sub_creating_time"] - tmp_cert["begin_time"]
        if feasible_sub.Status == 2:
            cut_expr, tmp_cert = insert_plane(model, feasible_sub, tmp_cert, feas=True,
                                              value_method=model.cbGetSolution)

        else:
            infeasible_sub.update()
            for variable in model._vars:
                infeasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                    round(model.cbGetSolution(variable))

            infeasible_sub.update()
            tmp_cert["infeas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"]
            infeasible_sub.optimize()
            if infeasible_sub.objVal == 0.:
                feasible_sub.computeIIS()
                feasible_sub.write("IIS.ilp")
                logging.debug("numerical problems in infeasible sub")
                exit()
            if infeasible_sub.Status != 2:
                infeasible_sub.write("infeasible_sub_infeasible.lp")
                logging.debug(infeasible_sub.Status)
                exit()

            tmp_cert["infeas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"] - tmp_cert["infeas_sub_creating_time"]

            cut_expr, tmp_cert = insert_plane(model, infeasible_sub, tmp_cert, feas=False,
                                              value_method=model.cbGetSolution)
            dump_if_over_time(model, certificate,
                              "Foo", False, target_directory, name_of_instance)
        expr_terms[iteration] = cut_expr
        model.cbLazy(cut_expr + model._epsilon <= 0)
        cb_subsolving.calls += 1
        cb_subsolving.tmp_calls += 1
cb_subsolving.calls = 0
cb_subsolving.tmp_calls = 0


def insert_plane(master, sub_problem, tmp_cert=dict(), feas = True,
                 value_method=lambda variable: variable.X):
    """
    Generate a certificate cut from a solved master and sub problem.

    Parameters
    ----------
    master : gurobipy.Model
        The master model (solved)
    sub_problem : gurobipy.Model
        The sub model (feasible or infeasible type) (solved)
    feas : bool, optional
        Is it a feasible or an infeasible sub?. The default is True.
        This means feasible.
    value_method : function, optional
        How variable values are evaluated.
        If used in callback, this has to be model.cbGetSolution.
        The default is lambda variable: variable.X.

    Returns
    -------
    cut_expr : gurobipy.LinExpr
        Containing the LHS of the certificate cut as gurobi expression.
    tmp_cert : dict
        Containing some other certificate information.

    """
    max_dual = 0.
    cut_expr = 0.
    tmp_cert["insert_plane_begin_time"] = time() - global_support.STARTING_TIME
    tmp_cert['feasible'] = feas
    tmp_cert['objective value'] = sub_problem.ObjVal
    logging.debug('normal vector:')
    tmp_cert['normal vector'] = dict()
    for variable in master._vars:
        coeff = sub_problem._constrs['fixing_' + variable.VarName].Pi
        max_dual = max(max_dual, abs(coeff))
        cut_expr += coeff*(variable - value_method(variable))
        if coeff != 0:
            tmp_cert['normal vector'].update({variable.VarName: coeff})
            logging.debug((variable.VarName, coeff))
    if max_dual >= 1e-06:
        cut_expr /= max_dual
    tmp_cert["norm"] = max_dual
    tmp_cert["insert_plane_time"] = time() - global_support.STARTING_TIME - tmp_cert["insert_plane_begin_time"]
    with open(os.path.join(global_support.TARGET_DIRECTORY,
                           global_support.NAME_OF_INSTANCE.split(".")[0] + "_intermediate.certificate"), 'rb+') as filehandle:
        filehandle.seek(-1, os.SEEK_END)
        filehandle.truncate()
    with open(os.path.join(global_support.TARGET_DIRECTORY,
                           global_support.NAME_OF_INSTANCE.split(".")[0] + "_intermediate.certificate"), 'a') as output_file:
        if tmp_cert["nbr"] == -1:
            output_file.write(str(f"\"{tmp_cert['nbr']}\":"))
        else:
            output_file.write(str(f",\"{tmp_cert['nbr']}\":"))
        output_file.write(json.dumps(tmp_cert))
        output_file.write("}")
    return cut_expr, tmp_cert


def dump_current_status(master, certificate, reduced_certificate, reduced):
    """Make a snapshot of the current status quo"""
    master.write('current_master.lp')
    certificate["reduced"] = reduced
    with open('current_certificate', 'w') as output_file:
        output_file.write(json.dumps(certificate))
    if reduced:
        reduced_certificate.write('reduced_certificate')


def initialize_problems(certificate, path_to_instance, target_directory, name_of_instance):
    """Read and solve the original problem. Copy it 3 times."""
    master = grb.read(path_to_instance)

    # optimize master problem firstly to determine optimal solution and first certificate point.
    master.setParam("MIPGap", 0)
    time_list = []
    time_list.append(time())
    time_before_solving = time()

    with open(os.path.join(global_support.TARGET_DIRECTORY,
                           global_support.NAME_OF_INSTANCE.split(".")[0] + "_intermediate.certificate"), 'a') as output_file:
              output_file.write(str("{}"))
    certificate["initial_solve"] = False
    certificate['certificate_without_reducing'] = False
    # solve the master problem for the first time.
    master.optimize(lambda model, where: cb_tracking(model, where, time_list))
    assert master.Status in [2, 3, 4, 5], "Master can not be optimized. Cancelling calculation."

    if master.Status == 2:
        for variable in master.getVars():
            variable._initial = variable.X
    # master.params.MIPGap = 0.9
    logging.info("Processed master problem successfully.")
    time_after_solving = time()
    certificate['time_to_solve'] = time_after_solving - time_before_solving
    certificate["initial_solve"] = True
    certificate['number_of_constraints'] = master.NumConstrs
    certificate['number_of_continuous_variables'] = master.NumVars - master.NumIntVars
    certificate['number_of_integral_variables'] = master.NumIntVars
    with open(os.path.join(target_directory,
                           name_of_instance.split(".")[0] + ".certificate"), 'w') as output_file:
        output_file.write(json.dumps(certificate))
    # master.setParam("logToConsole", 0)
    logging.debug("Time of optimal solution proof")
    logging.debug(- time_list[-1] + time())
    # master.setParam("MIPGap", 0)
    infeasible_sub = grb.read(path_to_instance)
    feasible_sub = grb.read(path_to_instance)
    infeasible_sub.setParam("QCPDual",1)
    feasible_sub.setParam("QCPDual",1)
    # master adaptation section
    # set output of gurobi models to "no output"
    master.setParam('LogToConsole', 0)
    feasible_sub.setParam('LogToConsole', 0)
    infeasible_sub.setParam('LogToConsole', 0)
    # safe the number of int variables as a number
    master._constrs = {}
    infeasible_sub._constrs = {}
    feasible_sub._constrs = {}
    # print logging message
    logging.info("problems initialized")
    return master, infeasible_sub, feasible_sub


def adjust_master(master, certificate):
    """
    Adjust the master problem.

    Remove all continuous variables and
    all constraints containing at least one continuous variable.
    Add epsilon variable.
    """
    number_of_kept_constraints = 0
    constraints_master = master.getConstrs()
    q_constraints_master = master.getQConstrs()
    gen_constraints_master = master.getGenConstrs()
    logging.debug('Creating master instance')
    remaining_master_constraints = len(constraints_master)
    for constr_nbr, constraint in enumerate(q_constraints_master + gen_constraints_master):
        master.remove(constraint)
    for constr_nbr, constraint in enumerate(constraints_master):
        leave_loop = False
        for variable in master.getVars():
            if (master.getCoeff(constraint, variable) != 0.
            and variable.VType == grb.GRB.CONTINUOUS):
                master.remove(constraint)
                remaining_master_constraints -= 1
                leave_loop = True
                break
        if not leave_loop:
            certificate['remaining_constraints']["remaining_constraint_{}".format(constr_nbr)] =\
                str(list(([master.getCoeff(constraint, var), var.varName]
                          for var in master.getVars()
                          if master.getCoeff(constraint, var) != 0.)))
            number_of_kept_constraints+=1

    master._vars = []
    number_of_kept_bounds = 0
    master._diameter = 0
    for variable in master.getVars():
        if variable.VType == grb.GRB.CONTINUOUS:
            master.remove(variable)
        else:
            master._vars.append(variable)
            if variable.lb != -grb.GRB.INFINITY:
                certificate["remaining_bounds"]["remaining_bound_{}".format(number_of_kept_bounds)] =\
                    variable.varName + " >= " + str(variable.lb)
                number_of_kept_bounds += 1
            if variable.ub != grb.GRB.INFINITY:
                certificate["remaining_bounds"]["remaining_bound_{}".format(number_of_kept_bounds)] =\
                    variable.varName + " <= " + str(variable.ub)
                number_of_kept_bounds += 1
            master._diameter = max(master._diameter, variable.ub - variable.lb)
    master.setObjective(0., grb.GRB.MINIMIZE)
    master._epsilon = master.addVar(lb=-1000, ub=1000.,
                                    name='epsilon', obj = -1.)
    logging.debug("Master Problem adaptations completed.")
    logging.debug('Creating infeasible_sub instance')
    certificate['number_of_kept_constraints'] = number_of_kept_constraints
    certificate['number_of_kept_bounds'] = number_of_kept_bounds


def adjust_inf_sub(infeasible_sub):
    """
    Adjust the infeasible-type sub problem.

    Add one variable for each constraint and bound (tmp) and a variable max_delta.
    Make the objective to min max_delta.
    Manipulate each constraint and bound that it can be violated by tmp.
    Add constraints tmp <= max_delta.
    Add dummy fixing constraints (new master solutions are inserted by
                                  manipulation of the RHS of the fixing
                                  constraints)
    """
    constraints_infeasible_sub = infeasible_sub.getConstrs()
    constraints_gen_infeasible_sub = infeasible_sub.getGenConstrs()

    constraints_q_infeasible_sub = infeasible_sub.getQConstrs()

    infeasible_sub.setObjective(0., grb.GRB.MINIMIZE)
    max_delta = infeasible_sub.addVar(name='max_delta', obj=1.)
    for constraint in constraints_infeasible_sub:
        if constraint.Sense == '<':
            tmp = infeasible_sub.addVar(column=grb.Column(-1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
        if constraint.Sense == '>':
            tmp = infeasible_sub.addVar(column=grb.Column(1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
        if constraint.Sense == '=':
            tmp = infeasible_sub.addVar(column=grb.Column(1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
            tmp = infeasible_sub.addVar(column=grb.Column(-1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)


    for constraint in constraints_gen_infeasible_sub:
        if constraint.Sense == '<':
            tmp = infeasible_sub.addVar(column=grb.Column(-1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
        if constraint.Sense == '>':
            tmp = infeasible_sub.addVar(column=grb.Column(1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
        if constraint.Sense == '=':
            tmp = infeasible_sub.addVar(column=grb.Column(1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)
            tmp = infeasible_sub.addVar(column=grb.Column(-1, constraint))
            infeasible_sub.addConstr(tmp <= max_delta)

    for constraint in list(constraints_q_infeasible_sub):
        if constraint.QCSense == '<':
            qLHS = infeasible_sub.getQCRow(constraint)
            qRHS = constraint.QCRHS
            qLHS.addTerms(-1., max_delta)
            infeasible_sub.addQConstr(qLHS<=qRHS)
            infeasible_sub.remove(constraint)
        if constraint.QCSense == '>':
            qLHS = infeasible_sub.getQCRow(constraint)
            qRHS = constraint.QCRHS
            qLHS.addTerms(1., max_delta)
            infeasible_sub.addQConstr(qLHS<=qRHS)
            infeasible_sub.remove(constraint)
        assert constraint.QCSense != '=', "modeling error quadractic equation"

    for variable in infeasible_sub.getVars():
        if variable.VType != grb.GRB.CONTINUOUS:
            infeasible_sub._constrs['fixing_' + variable.VarName] =\
                infeasible_sub.addConstr(variable == 0,
                                         name='fixing_' + variable.VarName)
            variable.VType = grb.GRB.CONTINUOUS
            continue
        if variable.LB != -grb.GRB.INFINITY:
            tmp = infeasible_sub.addVar()
            infeasible_sub.addConstr(tmp <= max_delta)
            infeasible_sub.addConstr(variable >= variable.LB - tmp)
            variable.LB = -grb.GRB.INFINITY
        if variable.UB != grb.GRB.INFINITY:
            tmp = infeasible_sub.addVar()
            infeasible_sub.addConstr(tmp <= max_delta)
            infeasible_sub.addConstr(variable <= variable.UB + tmp)
            variable.UB = grb.GRB.INFINITY
    logging.debug('Creating feasible_sub instance')


def adjust_feas_sub(feasible_sub, master, tmp_cert):
    """Adjust the feasible sub (see adjust_infeasible_sub) and solve it the first time."""
    for variable in feasible_sub.getVars():
        if variable.VType != grb.GRB.CONTINUOUS:
            feasible_sub._constrs['fixing_' + variable.VarName] =\
                feasible_sub.addConstr(variable == 0,
                                       name='fixing_' + variable.VarName)
            variable.VType = grb.GRB.CONTINUOUS
    logging.debug("Sub Problem adaptations completed. Going into loop...")
    logging.debug('Adding first cut')
    feasible_sub.update()
    logging.debug('Updating right hand side')
    tmp_cert["begin_time"] = time() - global_support.STARTING_TIME
    for variable in master._vars:
        feasible_sub._constrs['fixing_' + variable.VarName].RHS = variable.X
    feasible_sub.update()
    tmp_cert["feas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"]
    feasible_sub.optimize()
    tmp_cert["feas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"]


def hyperplane_computation(path_to_instance, name_of_instance,  target_directory="", args="",necessary_hyperplanes=0):
    """The main loop, see comments in code."""
    logging.debug(f'Path to instance: {path_to_instance}')
    open(os.path.join(target_directory,
                      name_of_instance.split(".")[0] + "_intermediate.certificate"),
         'w').close()

    certificate = dict()
    certificate['finished'] = "irregular"
    certificate['iterations'] = dict()
    certificate['remaining_constraints'] = dict()
    certificate['remaining_bounds'] = dict()
    certificate["reduced"] = False
    certificate['name_of_instance'] = name_of_instance
    with open(os.path.join(target_directory,
                               name_of_instance.split(".")[0] + ".certificate"), 'w') as output_file:
        output_file.write(json.dumps(certificate))
    starting_time = time()
    global_support.STARTING_TIME = starting_time

    master, infeasible_sub, feasible_sub = initialize_problems(certificate,
                                                               path_to_instance, target_directory, name_of_instance)

    if master.Status in [3, 4, 5]:
        logging.info("Problem is infeasible or unbounded. Returning empty certificate.")
        certificate['feasibility'] = False
        certificate['terminated_correctly'] = True
        dump_if_over_time(master, certificate,
                          "Foo", False,target_directory, name_of_instance)
        return certificate
    else:
        logging.info("Problem is feasible. Proceding normally.")
        certificate['feasibility'] = True
        dump_if_over_time(master, certificate,
                          "Foo", False,target_directory, name_of_instance)
    adjust_master(master, certificate)
    with open(os.path.join(target_directory,
                           name_of_instance.split(".")[0] + ".certificate"), 'w') as output_file:
        output_file.write(json.dumps(certificate))
    logging.info("master adjusted")

    tmp_cert = {"nbr": -1}

    adjust_inf_sub(infeasible_sub)
    logging.info("sub inf adjusted")
    adjust_feas_sub(feasible_sub, master, tmp_cert)
    logging.info("sub feas adjusted")

    certificate_calculation_begin_time = time()
    tmp_cert['suspension point'] = {variable.VarName: variable.X
                                            for variable in master._vars
                                            if variable.X != 0}
    cut_expr, tmp_cert = insert_plane(master, feasible_sub, tmp_cert)
    master.addConstr(cut_expr + master._epsilon <= 0, name="iteration_"+str(-1))
    necessary_hyperplanes+=1
    certificate['iterations'][f"iteration{-1}"] = tmp_cert.copy()

    dump_if_over_time(master, certificate,
                      "Foo", False,target_directory, name_of_instance)

    if global_support.CALLBACK:
        master.setParam("LazyConstraints", 1)
    else:
        master.setParam("MIPGap", global_support.OPT_MODE)

    expr_terms = {}
    logging.debug('Entering the first main loop')
    subs_calls_old = 0
    if args.boxing and master._diameter not in [0., 1.]:
        if master._diameter == float("inf"):
            my_gen = [10**i for i in range(global_support.BOXING_PARAM)]
        else:

            my_gen = list(set([int(int(master._diameter) ** (float(i)/(global_support.BOXING_PARAM - 1))) for i in range(global_support.BOXING_PARAM - 1)]))
        for distance in my_gen:
            box_loop(target_directory, name_of_instance,
                     master, feasible_sub, infeasible_sub, expr_terms,
                     certificate)
    tmp_obj_value = 1000.
    while tmp_obj_value != 0.:
        subs_calls_old = cb_subsolving.calls
        if global_support.CALLBACK:
            master.optimize(lambda model, where: cb_subsolving(model, where,
                                                               feasible_sub,
                                                               infeasible_sub,
                                                               certificate,
                                                               expr_terms,
                                                               target_directory,
                                                               name_of_instance))
            dump_if_over_time(master, certificate,
                              "Foo", False,target_directory, name_of_instance)
        else:
            master.optimize()
            logging.info("Iteration: " + str(cb_subsolving.calls) + ", Current Progress: "
                         + str(- master.objval))
            iteration = cb_subsolving.calls
            tmp_cert = certificate['iterations']['iteration' + str(iteration)] = {"nbr": iteration}
            tmp_cert["begin_time"] = time() - global_support.STARTING_TIME
            tmp_cert['suspension point'] = {variable.VarName: variable.X
                                                for variable in master._vars
                                                if variable.X != 0}
            feasible_sub.update()
            dump_if_over_time(master, certificate,
                              "Foo", False,target_directory, name_of_instance)
            for variable in master._vars:
                feasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                    variable.X
            feasible_sub.update()
            tmp_cert["feas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"]
            feasible_sub.optimize()
            tmp_cert["feas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["feas_sub_creating_time"] - tmp_cert["begin_time"]
            if feasible_sub.Status == 2:
                cut_expr, tmp_cert = insert_plane(master, feasible_sub, tmp_cert, feas=True)

            else:
                infeasible_sub.update()
                for variable in master._vars:
                    infeasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                        variable.X

                infeasible_sub.update()
                tmp_cert["infeas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"]
                infeasible_sub.optimize()
                tmp_cert["infeas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"] - tmp_cert["infeas_sub_creating_time"]

                cut_expr, tmp_cert = insert_plane(master, infeasible_sub, tmp_cert, feas=False)
            expr_terms[iteration] = cut_expr
            cb_subsolving.calls += 1
            cb_subsolving.tmp_calls += 1

        if master.status == 2 and - master.objval <= 1e-06:
            exit_cond = True
        else:
            exit_cond = False
        logging.debug("leaving main loop")
        for i in range(subs_calls_old, cb_subsolving.calls):
            master.addConstr(expr_terms[i] + master._epsilon <= 0, name="iteration_"+str(i))
        master.update()
        if exit_cond:
            break
    master.write(os.path.join(target_directory, name_of_instance + '_certificate' + '.lp'))

    certificate_calculation_end_time = time()
    certificate["calculation_time"] = (certificate_calculation_end_time -
                                       certificate_calculation_begin_time)

    certificate['number_of_necessary_hyperplanes'] = necessary_hyperplanes
    certificate['reduce_time'] = -1.
    certificate['certificate_without_reducing'] = True
    dump_if_over_time(master, certificate,
                      "Foo", False,target_directory, name_of_instance)
    with open(os.path.join(target_directory,
                           name_of_instance.split(".")[0] + "_before_reducing.certificate"), 'w') as output_file:
        output_file.write(json.dumps(certificate))
    with open(os.path.join(target_directory,
                        name_of_instance.split(".")[0] + ".certificate"), 'w') as output_file:
        output_file.write(json.dumps(certificate))
    if args.reduce:
        reduce_certificate(master, certificate, target_directory,
                           name_of_instance)
    cb_subsolving.calls = 0
    return certificate


def box_loop(target_directory, name_of_instance,master, feasible_sub,
             infeasible_sub, expr_terms, certificate,
             distance=1.):
    """Create constraints around optimal solutions of master. Call one optimize with callback."""
    logging.debug('Creating box')
    box_constraints = []
    for variable in master._vars:
        box_constraints.append(master.addConstr(variable >= variable._initial - distance,
                                                name='box_lower_' + variable.VarName))
        box_constraints.append(master.addConstr(variable <= variable._initial + distance,
                                                name='box_upper_' + variable.VarName))
    master.update()
    while True:
        cb_subsolving_calls_old = cb_subsolving.calls
        if global_support.CALLBACK:
            dump_if_over_time(master, certificate,
                              "Foo", False,target_directory, name_of_instance)
            master.optimize(lambda model, where: cb_subsolving(model, where,
                                                               feasible_sub,
                                                               infeasible_sub,
                                                               certificate,
                                                               expr_terms,
                                                               target_directory,
                                                               name_of_instance))
        else:
            master.optimize()
            dump_if_over_time(master, certificate,
                              "Foo", False, target_directory, name_of_instance)
            logging.info("Iteration: " + str(cb_subsolving.calls) + ", Current Progress: "
                         + str(- master.objval))
            iteration = cb_subsolving.calls
            tmp_cert = certificate['iterations']['iteration' + str(iteration)] = {"nbr": iteration}
            tmp_cert["begin_time"] = time() - global_support.STARTING_TIME
            tmp_cert['suspension point'] = {variable.VarName: variable.X
                                                for variable in master._vars
                                                if variable.X != 0}
            feasible_sub.update()
            for variable in master._vars:
                feasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                    variable.X
            feasible_sub.update()
            tmp_cert["feas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"]
            feasible_sub.optimize()
            tmp_cert["feas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["feas_sub_creating_time"] - tmp_cert["begin_time"]
            if feasible_sub.Status == 2:
                cut_expr, tmp_cert = insert_plane(master, feasible_sub, tmp_cert, feas=True)

            else:
                infeasible_sub.update()
                for variable in master._vars:
                    infeasible_sub._constrs['fixing_' + variable.VarName].RHS =\
                        variable.X

                infeasible_sub.update()
                tmp_cert["infeas_sub_creating_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"]
                infeasible_sub.optimize()
                tmp_cert["infeas_sub_opt_time"] = time() - global_support.STARTING_TIME - tmp_cert["begin_time"] - tmp_cert["feas_sub_creating_time"] - tmp_cert["feas_sub_opt_time"] - tmp_cert["infeas_sub_creating_time"]

                cut_expr, tmp_cert = insert_plane(master, infeasible_sub, tmp_cert, feas=False)
            expr_terms[iteration] = cut_expr
            cb_subsolving.calls += 1
            cb_subsolving.tmp_calls += 1
        if master.status == 2 and - master.objval <= 1e-06:
            exit_cond = True
        else:
            exit_cond = False

        logging.debug("leaving box loop")
        for i in range(cb_subsolving_calls_old, cb_subsolving.calls):
            master.addConstr(expr_terms[i] + master._epsilon <= 0, name="iteration_"+str(i))
        master.update()
        if exit_cond:
            break
    logging.debug('Removing artificial box')
    for hyperplane in box_constraints:
        master.remove(hyperplane)
    master.update()


def reduce_certificate(master, certificate, target_directory, name_of_instance):
    """Reduce the certificate by successively turning constraints."""
    reduce_begin = time()
    certificate_size = 0
    removed_constraints = 0
    lp_removed_constraints = 0
    uncertain_keeps = 0
    logging.debug('removing')
    logging.debug("iterations: " + str(cb_subsolving.calls))
    logging.debug('removing')
    necessary_hyperplanes = 0.
    for hyperplane in [hp for hp in master.getConstrs()
                       if hp.ConstrName.startswith("iteration")]:
        logging.debug(f"Time_usage_currently_{time()-global_support.STARTING_TIME}")
        dump_if_over_time(master, certificate,
                          "foo", False,
                          target_directory,
                          name_of_instance)

        certificate_size += 1
        hyperplane.Sense = ">"
        master.chgCoeff(hyperplane, master._epsilon, 0.)
        time_list = []
        time_list.append(time())
        master.optimize(lambda model, where: cb_tracking(model, where, time_list))

        if master.Status != grb.GRB.OPTIMAL:
            tmp_bool = True
        else:
            tmp_bool = (master.ObjVal >= -1e-09)
        if tmp_bool:
            try:
                del certificate['iterations'][hyperplane.ConstrName.replace('_', '')]
            except KeyError:
                logging.debug("Removing hyperplane failed.")
            master.remove(hyperplane)
            removed_constraints += 1
            logging.info("removed hyperplane " +
                         hyperplane.ConstrName +
                         " after " +
                         str(master.NodeCount) +
                         " MIP Nodes")
            if master.NodeCount == 0.0:
                lp_removed_constraints += 1
        else:
            hyperplane.Sense = "<"
            master.chgCoeff(hyperplane, master._epsilon, 1.)
            logging.info("kept " + hyperplane.ConstrName +
                         " with optimal value " + str(master.ObjVal))
            necessary_hyperplanes+=1
            if -1e-04 <= master.ObjVal:
                uncertain_keeps += 1
    certificate['number_of_necessary_hyperplanes'] = necessary_hyperplanes
    master.write(os.path.join(target_directory,
                              name_of_instance + '_reduced_certificate' + '.lp'))
    msg_str = f"Removed {lp_removed_constraints} (LP) "
    msg_str += f"{removed_constraints} (Total) of "
    certificate["nbr_reduced_constraints"] = removed_constraints
    msg_str += f"{certificate_size} certificate points. "
    certificate["nbr_planes_before_reducing"] = certificate_size
    msg_str += f"{uncertain_keeps} of remaining points uncertain."
    logging.info(msg_str)
    reduce_end = time()
    logging.info(time_list)
    certificate["reduce_time"] = (reduce_end - reduce_begin)


def main():
    parser = argparse.ArgumentParser(description='Compute a certificate.')
    parser.add_argument("instance")
    parser.add_argument(
        '-d', '--debug',
        help="Print lots of debugging statements",
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        help="Be verbose",
        action="store_const", dest="loglevel", const=logging.INFO,
    )
    parser.add_argument(
        '-l', '--local',
        action="store_true",
        help='Is this file a local lp-file. Crashes with -p'
    )
    parser.add_argument(
        '-p', '--mps',
        help="Computes local mps file. Crashes with -l",
        action="store_true",
    )
    parser.add_argument(
        '-r', '--reduce',
        action="store_true",
        help='Do you want to reduce the certificate?'
    )
    parser.add_argument(
        '-b', '--boxing',
        action="store_true",
        help='Do you want to search for certificate points in box around optimal solution first?'
    )
    parser.add_argument(
        '-m', '--multiple',
        help="Activate to type a directory name and calculate"+\
            " a certificate for all instance files in folder.",
        action="store_true",
    )
    parser.add_argument(
        '-f', '--file_path',
        help="give path to file",
        action="store_true",
    )
    parser.add_argument(
        '-c', '--compute_dir',
        help="compute in compute_dir",
        action="store_true",
    )
    parser.add_argument(
        '--opt_mode',
        help="If == 1, then optimize master with callbacks. ElIf \in [0,1), optimize with MIPGap args.opt_mode. Else Error.",
        default=1.,
        type=float,
    )
    parser.add_argument(
        '--memory_limit',
        help="Memory limit in MiB. Default 5MiB.",
        default=5000.,
        type=int,
    )
    parser.add_argument(
        '--boxing_param',
        help="How many boxing loops. Box size grows exponentially, based on diameter max(ub-lb) of integer variables of problem. Default 4.",
        default=4,
        type=int,
    )
    parser.add_argument(
        '--reset_param',
        help="Restarts Callback Loop after --reset_param iterations. Default infinity.",
        default=float("inf"),
        type=float,
    )
    parser.add_argument(
        '--time_limit',
        help="Time limit in seconds. Default 3500.",
        default=3500,
        type=int,
    )
    parser.add_argument(
        '--hpc',
        help="Less data produced",
        type=bool,
    )
    args = parser.parse_args()
    arg_str = "_"
    arg_str += "boxing_" + str(args.boxing) + "_"
    arg_str += "opt_mode_" + str(args.opt_mode) + "_"
    arg_str += "boxing_param_" + str(args.boxing_param) + "_"
    arg_str += "reset_param_" + str(args.reset_param) + "_"
    arg_str += "time_limit_" + str(args.time_limit)

    global_support.MEMORY_LIMIT = float(args.memory_limit)
    global_support.TIME_LIMIT = float(args.time_limit)
    global_support.RESET = float(args.reset_param)
    global_support.BOXING_PARAM = int(args.boxing_param)
    if float(args.opt_mode) != 1.:
        global_support.CALLBACK = False
    global_support.OPT_MODE = args.opt_mode

    if not os.path.exists("Loggingfiles"):
        os.mkdir("Loggingfiles")
    logname = os.path.join("Loggingfiles", f"{args.instance}.log")

    logging.basicConfig(level=args.loglevel,  filemode="w", format='%(asctime)s, %(levelname)s: %(message)s', stream=sys.stdout)

    instance_list = [args.instance.split('.')[0]]

    if args.local and not args.multiple:
        path_to_instance_list = [os.path.join(os.getcwd(), str(instance_list[0]) +'.lp')]
    if args.mps and not args.multiple:
        path_to_instance_list = [os.path.join(os.getcwd(), str(instance_list[0]) +'.mps')]
    if args.file_path and not args.multiple:
        logging.debug(f'instance_list:{instance_list}')
        instance_list = [os.path.basename(args.instance).split('.')[0]]
        logging.debug(f'instance_list:{instance_list}')
        path_to_instance_list = [args.instance]
    elif not args.local and not args.multiple:
        downloader(instance_list[0])
        path_to_instance_list = [os.path.join(os.getcwd(), 'MIPLIBing_cache',
                                              'MIPLIB2017_Collection', instance_list[0] + '.mps')]
    elif args.multiple:
        path_to_instance_list = list((os.path.join(os.getcwd(), args.instance, file_name)
                                      for file_name in os.listdir(os.path.join(os.getcwd(),
                                                                               args.instance))))
        instance_list = list((file_name.split(".")[0]
                              for file_name in os.listdir(os.path.join(os.getcwd(),
                                                                       args.instance))))

    for instance, path_to_instance in zip(instance_list, path_to_instance_list):
        global_support.NAME_OF_INSTANCE = instance
        target_directory = os.path.join(os.getcwd(), 'Certificates', instance)
        if args.hpc:
            hpc_name = args.instance.split('/')[-2]
            hpc_name += arg_str
            target_directory = os.path.join(os.getcwd(), \
                    'computations_'+ now.strftime("%y_%m_%d") + '_' + hpc_name)
        global_support.TARGET_DIRECTORY = target_directory
        if not os.path.exists(target_directory):
            os.mkdir(target_directory)
        if args.hpc:
            pass
        else:
            copy(path_to_instance, target_directory)
        logging.debug(f'instance:{instance}')
        logging.debug(f'target_directory:{target_directory}')

        certificates = hyperplane_computation(path_to_instance,
                                              instance,
                                              target_directory = target_directory,
                                              args=args)

        certificates["finished"] = "regular"
        with open(os.path.join(target_directory,
                               instance.split(".")[0] + ".certificate"), 'w') as output_file:
            output_file.write(json.dumps(certificates))

    if args.loglevel <= 20:
        open("todelete.cprofile", "w").close()
        f = open(os.path.join("Loggingfiles", f"{args.instance}.cprofile"), 'w')
        sys.stdout = f
    else:
        f = open("todelete.cprofile", "w")
        sys.stdout = f



if __name__ == "__main__":
    main()
