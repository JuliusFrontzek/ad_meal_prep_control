import do_mpc


class StateEstimator(do_mpc.estimator.Estimator):
    """
    For now this functions as a state-feedback but will facilitate the EKF
    as well as add some noise to the measured variables.
    """

    def __init__(self, model):
        super().__init__(model)

    def make_step(self, y0):
        """Return the measurement ``y0``."""
        return y0
