import cdutil
import cdp.cdp_metric

class Averager(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(Averager, self).__init__(metric_path)

    def compute(self, variable, axes, weights):
        return cdutil.averager(variable, axis=axes, weights=weights)