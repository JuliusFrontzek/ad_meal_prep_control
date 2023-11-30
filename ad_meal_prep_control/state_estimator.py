import do_mpc
import numpy as np
from pathlib import Path


class StateFeedback(do_mpc.estimator.Estimator):
    """
    For now this functions as a state-feedback but will facilitate the EKF
    as well as add some noise to the measured variables.
    """

    def __init__(self, model, simulator):
        super().__init__(model)
        self._simulator = simulator

    def make_step(self, y) -> np.ndarray:
        """Return the actual state of the system. The parameter y is passed as a dummy but not actually used to match the signature of the MHE 'make_step' method."""
        return np.array(self._simulator.x0.master)


def mhe_setup(
    model: do_mpc.model.Model,
    t_step: float,
    n_horizon: int,
    x_ch_in: np.ndarray,
    x_pr_in: np.ndarray,
    x_li_in: np.ndarray,
    P_x: np.ndarray,
    P_v: np.ndarray,
    ch4_outflow_rate: np.ndarray,
    store_full_solution: bool,
    hsllib: Path = None,
) -> do_mpc.estimator.MHE:
    """
    MHE setup
    """
    num_states = model._x.size
    mhe = do_mpc.estimator.MHE(
        model, p_est_list=[f"x_{i+1}" for i in range(num_states)]
    )
    setup_mhe = {
        "t_step": t_step,
        "n_horizon": n_horizon,
        "store_full_solution": store_full_solution,
        "meas_from_data": True,
    }

    if hsllib is not None:
        setup_mhe["nlpsol_opts"] = {
            "ipopt.linear_solver": "MA27",
            "ipopt.hsllib": str(hsllib),
        }

    mhe.set_param(**setup_mhe)

    mhe.set_default_objective(P_x=P_x, P_v=P_v)

    if num_states == 20:
        tvp_template = mhe.get_tvp_template()

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            n_steps = min(t_now_idx, n_horizon)

            for k in range(n_steps):
                tvp_template["_tvp", k, "v_ch4_dot_out"] = ch4_outflow_rate[
                    t_now_idx + k
                ]

            return tvp_template

        mhe.set_tvp_fun(tvp_fun)

    p_num = mhe.get_p_template()
    p_num["x_ch_in"] = np.array(
        [np.mean(x_ch_in[:, i]) for i in range(x_ch_in.shape[1])]
    )
    p_num["x_pr_in"] = np.array(
        [np.mean(x_pr_in[:, i]) for i in range(x_pr_in.shape[1])]
    )
    p_num["x_li_in"] = np.array(
        [np.mean(x_li_in[:, i]) for i in range(x_li_in.shape[1])]
    )

    def p_fun(t_now):
        return p_num

    mhe.set_p_fun(p_fun)

    # Every state but the 13th state must not be negative
    for i in range(1, num_states + 1):
        if i != 13:
            mhe.bounds["lower", "_x", f"x_{i}"] = 0.0

    mhe.setup()

    return mhe
