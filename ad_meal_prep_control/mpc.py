import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from ad_meal_prep_control.params_R3 import *
from ad_meal_prep_control.utils import (
    CostFunction,
    Bound,
    NlConstraint,
    LimitedSubstrate,
    SetpointFunction,
    Disturbances,
)
from pathlib import Path

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def mpc_setup(
    *,
    model: do_mpc.model.Model,
    n_horizon: int,
    n_robust: int,
    t_step: float,
    xi_ch_norm: np.ndarray,
    xi_pr_norm: np.ndarray,
    xi_li_norm: np.ndarray,
    compile_nlp: bool,
    ch4_outflow_rate: np.ndarray,
    cost_func: CostFunction,
    substrate_costs: list[float],
    substrate_cost_formulation: str,
    store_full_solution: bool,
    disturbances: Disturbances,
    gas_storage_bound_fraction: float,
    bounds: list[Bound] = None,
    nl_cons: list[NlConstraint] = None,
    rterm: str = None,
    hsllib: Path = None,
    suppress_ipopt_output: bool = False,
    ch4_set_point_function: SetpointFunction = None,
    theta: np.ndarray,
    limited_substrates: list[LimitedSubstrate] = None,
) -> do_mpc.controller.MPC:
    num_states = model._x.size

    mpc = do_mpc.controller.MPC(model)

    setup_mpc = {
        "n_horizon": n_horizon,
        "n_robust": n_robust,
        "open_loop": False,
        "t_step": t_step,
        "state_discretization": "collocation",
        "collocation_type": "radau",
        "collocation_deg": 2,
        "collocation_ni": 1,
        "store_full_solution": store_full_solution,
        # "nl_cons_check_colloc_points": True,
    }

    if hsllib is not None:
        # Use MA27 linear solver in ipopt for faster calculations:
        setup_mpc["nlpsol_opts"] = {
            "ipopt.linear_solver": "MA27",
            "ipopt.hsllib": str(hsllib),
        }
    else:
        setup_mpc["nlpsol_opts"] = {}

    # setup_mpc["nlpsol_opts"] = {
    #     "ipopt.linear_solver": "spral",
    #     "ipopt.hsllib": str(hsllib),
    # }

    if suppress_ipopt_output:
        setup_mpc["nlpsol_opts"]["ipopt.print_level"] = 0

    mpc.set_param(**setup_mpc)

    tvp_template = mpc.get_tvp_template()

    def dictated_sub_tvp_setup(t_now: float):
        if disturbances.dictated_feeding is not None:
            for feed_idx, dictated_sub_properties in enumerate(
                disturbances.dictated_feeding.values()
            ):
                for k in range(n_horizon + 1):
                    for dictated_sub in dictated_sub_properties:
                        start_time = dictated_sub[0]
                        end_time = dictated_sub[1]
                        feed_amount = dictated_sub[2]

                        time_k = t_now + t_step * k
                        if time_k >= start_time and time_k < end_time:
                            tvp_template[
                                "_tvp", k, "dictated_sub_feed", feed_idx
                            ] = feed_amount
                            break
                        else:
                            tvp_template["_tvp", k, "dictated_sub_feed", feed_idx] = 0.0

    if num_states == 20:  # i.e. if we consider the gas storage
        mpc.set_nl_cons(
            "max_vol_gas_storage",
            model.x["x_19"] + model.x["x_20"] + model.aux["V_H2O"] / V_GAS_STORAGE_MAX,
            ub=1.0 - gas_storage_bound_fraction,
            soft_constraint=True,
            penalty_term_cons=1e1,
            maximum_violation=gas_storage_bound_fraction,
        )

        mpc.set_nl_cons(
            "min_vol_gas_storage",
            -(
                model.x["x_19"]
                + model.x["x_20"]
                + model.aux["V_H2O"] / V_GAS_STORAGE_MAX
            ),
            ub=-gas_storage_bound_fraction,
            soft_constraint=True,
            penalty_term_cons=1e1,
            maximum_violation=gas_storage_bound_fraction,
        )

        # for i in range(2):
        #     mpc.set_nl_cons(
        #         f"x_{19+i}_lower_bound",
        #         -model.x[f"x_{19+i}"],
        #         ub=-0.05,
        #         soft_constraint=True,
        #         penalty_term_cons=1e4,
        #         maximum_violation=0.05,
        #     )
        mpc.bounds["lower", "_x", "x_19"] = 0.0
        mpc.bounds["lower", "_x", "x_20"] = 0.0

        # mpc.set_nl_cons(
        #     f"pH_lower_bound",
        #     -model.aux[f"y_4"],
        #     ub=-6.0,
        #     soft_constraint=True,
        #     penalty_term_cons=1e3,
        #     maximum_violation=0.2,
        # )

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            mean_ch4_outflow_rate = np.mean(
                ch4_outflow_rate
            )  # [t_now_idx : t_now_idx + n_horizon]
            for k in range(n_horizon + 1):
                tvp_template["_tvp", k, "v_ch4_dot_tank_out"] = ch4_outflow_rate[
                    t_now_idx + k
                ]
                tvp_template[
                    "_tvp", k, "v_ch4_dot_tank_out_mean"
                ] = mean_ch4_outflow_rate
                tvp_template["_tvp", k, "theta"] = theta

            if ch4_set_point_function is not None:
                tvp_template[
                    "_tvp", :, "v_ch4_dot_tank_in_setpoint"
                ] = ch4_set_point_function.get_current_setpoint(t_now)

            dictated_sub_tvp_setup(t_now)

            return tvp_template

    else:

        def tvp_fun(t_now):
            for k in range(n_horizon + 1):
                tvp_template["_tvp", k, "theta"] = theta
            if ch4_set_point_function is not None:
                tvp_template[
                    "_tvp", :, "v_ch4_dot_tank_in_setpoint"
                ] = ch4_set_point_function.get_current_setpoint(t_now)

            dictated_sub_tvp_setup(t_now)

            return tvp_template

    mpc.set_tvp_fun(tvp_fun)

    mpc.set_objective(lterm=eval(cost_func.lterm), mterm=eval(cost_func.mterm))

    if substrate_cost_formulation in ["linear", "quadratic"]:
        substrate_costs = np.array(substrate_costs)
        substrate_costs /= np.max(substrate_costs[substrate_costs > 0])

        sub_cost_rterms = []

        if substrate_cost_formulation == "linear":
            for idx, cost in enumerate(substrate_costs):
                sub_cost_rterms.append(f"{cost} * (model.u['u_norm'][{idx}])")
        else:
            for idx, cost in enumerate(substrate_costs):
                sub_cost_rterms.append(f"{cost} * (model.u['u_norm'][{idx}])**2")

        sub_cost_rterm = " + ".join(sub_cost_rterms)

        if rterm is None:
            rterm = sub_cost_rterm
        else:
            rterm += f" + {sub_cost_rterm}"

    if rterm is not None:
        mpc.set_rterm(rterm=eval(rterm))

    # Hard constraints
    mpc.bounds["lower", "_u", "u_norm"] = 0.0
    mpc.bounds["upper", "_u", "u_norm"] = 1.0

    if bounds is not None:
        for bound in bounds:
            mpc.bounds[
                bound.direction, bound.variable_type, bound.variable
            ] = bound.value

    # # Soft constraints with slack variables
    if nl_cons is not None:
        for idx, nl_con in enumerate(nl_cons):
            mpc.set_nl_cons(
                expr_name=f"nl_cons_{idx}",
                expr=eval(nl_con.expression),
                ub=nl_con.ub,
                soft_constraint=nl_con.soft_constraint,
                penalty_term_cons=nl_con.penalty_term_cons,
                maximum_violation=nl_con.maximum_violation,
            )

    mpc.set_uncertainty_values(
        xi_ch_norm=xi_ch_norm, xi_pr_norm=xi_pr_norm, xi_li_norm=xi_li_norm
    )

    mpc.setup()
    if compile_nlp:
        mpc.compile_nlp(
            overwrite=True, compiler_command="gcc -fPIC -shared -O3 nlp.c -o nlp.so"
        )

    return mpc
