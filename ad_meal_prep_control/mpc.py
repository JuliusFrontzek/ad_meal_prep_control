import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from params_R3 import *
from ad_meal_prep_control.utils import CostFunction, Bound, NlConstraint
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
    x_ch_in: np.ndarray,
    x_pr_in: np.ndarray,
    x_li_in: np.ndarray,
    compile_nlp: bool,
    vol_flow_rate: np.ndarray,
    cost_func: CostFunction,
    substrate_costs: list[float],
    consider_substrate_costs: bool,
    bounds: list[Bound] = None,
    nl_cons: list[NlConstraint] = None,
    rterm: str = None,
    hsllib: Path = None,
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
    }

    if hsllib is not None:
        # Use MA27 linear solver in ipopt for faster calculations:
        setup_mpc["nlpsol_opts"] = {
            "ipopt.linear_solver": "MA27",
            "ipopt.hsllib": str(hsllib),  # "./coinhsl*/builddir/libcoinhsl.so",
        }

    mpc.set_param(**setup_mpc)

    mpc.set_param(store_full_solution=True)

    if num_states == 20:  # i.e. if we consider the gas storage
        tvp_template = mpc.get_tvp_template()

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            for k in range(n_horizon + 1):
                tvp_template["_tvp", k, "v_ch4_dot_out"] = vol_flow_rate[t_now_idx + k]

            return tvp_template

        mpc.set_nl_cons(
            "max_vol_gas_storage",
            model._aux_expression["v_gas_storage"],
            ub=V_GAS_STORAGE_MAX,
            soft_constraint=True,
            penalty_term_cons=1e7,
        )
        mpc.bounds["lower", "_x", "x_19"] = 0.0
        mpc.bounds["lower", "_x", "x_20"] = 0.0

        mpc.set_tvp_fun(tvp_fun)

    mpc.set_objective(lterm=eval(cost_func.lterm), mterm=eval(cost_func.mterm))

    if consider_substrate_costs:
        substrate_costs = np.array(substrate_costs)
        substrate_costs /= np.min(substrate_costs[substrate_costs > 0])

        sub_cost_rterms = []
        for idx, cost in enumerate(substrate_costs):
            sub_cost_rterms.append(f"{cost} * model.u['u_norm'][{idx}]**2")

        sub_cost_rterm = " + ".join(sub_cost_rterms)

        if rterm is None:
            rterm = sub_cost_rterm
        else:
            rterm += sub_cost_rterm

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

    mpc.set_uncertainty_values(x_ch_in=x_ch_in, x_pr_in=x_pr_in, x_li_in=x_li_in)

    mpc.setup()
    if compile_nlp:
        mpc.compile_nlp(
            overwrite=True, compiler_command="gcc -fPIC -shared -O3 nlp.c -o nlp.so"
        )

    return mpc
