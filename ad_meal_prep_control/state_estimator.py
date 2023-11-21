import do_mpc


class StateEstimator(do_mpc.estimator.Estimator):
    """
    For now this functions as a state-feedback but will facilitate the EKF
    as well as add some noise to the measured variables.
    """

    def __init__(self, model):
        super().__init__(model)

    def estimate_x(self, y):
        """Return the measurement ``y`` because this is a state estimator."""
        return y
