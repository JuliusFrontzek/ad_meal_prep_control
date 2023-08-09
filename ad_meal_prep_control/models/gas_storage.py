from casadi import *
from casadi.tools import *
import sys
import os
from dataclasses import dataclass

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


@dataclass
class GasConditions:
    p0: float
    p_gas_storage: float
    p_norm: float
    T_norm: float
    p_co2_phase_change: float
    p_ch4_phase_change: float
    p_h2o: float
    T: float


def gas_storage_model(cons: GasConditions):
    """
    Gas storage model with constant vapour pressure of water.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")
    num_inputs = 2  # V´_g and V´_CH4_out

    # Input
    u = model.set_variable(var_type="_u", var_name="u", shape=(num_inputs, 1))

    # States
    x = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(2)
    ]

    v_h2o = model.set_expression(
        "V_H2O",
        1.0
        / (1.0 - cons.p_h2o / cons.p_gas_storage)
        * cons.p_h2o
        / cons.p_gas_storage
        * (x[0] + x[1]),
    )

    y_co2 = model.set_expression(
        "y_co2",
        x[1] / (x[0] + x[1] + v_h2o),
    )

    y_h2o = model.set_expression(
        "y_h2o",
        v_h2o / (x[0] + x[1] + v_h2o),
    )

    # Differential equations
    model.set_rhs(
        "x_1",
        u[0]
        * cons.p_norm
        / cons.p_gas_storage
        * cons.T
        / cons.T_norm
        * cons.p_ch4_phase_change
        / (cons.p_ch4_phase_change + cons.p_co2_phase_change + cons.p_h2o)
        - u[1],
    )  # V_CH4

    model.set_rhs(
        "x_2",
        u[0]
        * cons.p_norm
        / cons.p_gas_storage
        * cons.T
        / cons.T_norm
        * cons.p_co2_phase_change
        / (cons.p_ch4_phase_change + cons.p_co2_phase_change + cons.p_h2o)
        - y_co2
        / (1.0 - y_co2)
        / (1.0 - y_co2 * y_h2o / ((1.0 - y_co2) * (1.0 - y_h2o)))
        * (u[1] * (1.0 + y_h2o / (1.0 - y_h2o))),
    )  # V_CO2

    # Build the model
    model.setup()

    return model
