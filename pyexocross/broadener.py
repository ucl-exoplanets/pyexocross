
class Broadener:

    def __init__(self, ratio=1.0):

        self._ratio = ratio


    def set_temperature_pressure(self, T, P):
        self.T = T
        self.P = P



    @property
    def species(self):
        raise NotImplementedError

    @property
    def ratio(self):
        return self._ratio
    
    @ratio.setter
    def ratio(self, value):
        self._ratio = value

    def calculate_gamma(self, transitions):
        raise NotImplementedError


    def gamma(self, transitions):
        return self.calculate_gamma(transitions)*self.ratio
 