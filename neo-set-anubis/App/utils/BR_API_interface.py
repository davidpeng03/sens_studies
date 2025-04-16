import sys,os
sys.path.append(os.getcwd())

from Pipeline.BR_calculator import BRInterface


class BR_API_Interface:
    
    def __init__(self, x_min, y_min, x_max, y_max, step, x_vals, calculation_type, particles=None, model="HNL", params={"V1": 1, "V2": 1, "V3": 1}, masses={"N1": 1}):
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max
        self.step = step
        self.x_vals = x_vals
        self.channels = particles if particles is not None else []  # List of channels
        self.calculation_type = calculation_type
        self.br_calculator = BRInterface()
        self.br_calculator.set_calculation_method("Python", True)
        self.br_calculator.set_model(model)
        self.br_calculator.set_params(params)
        self.br_calculator.set_masses(masses)
        self.params = params
        self.masses = masses
        self.y_vals = {tuple(channel): [] for channel in self.channels}
        self._calculate_y_vals()

    def _calculate_y_vals(self):
        for channel in self.channels:
            y_vals_channel = []
            for x in self.x_vals:
                self.br_calculator.set_masses({"N1": x})
                y_vals_channel.append(self.br_calculator.calculate(self.calculation_type, "N1", channel))
            self.y_vals[tuple(channel)] = y_vals_channel

    def set_x_min(self, x_min):
        if x_min > self.x_min:
            self.x_vals = [x for x in self.x_vals if x >= x_min]
        else:
            new_vals = list(range(x_min, self.x_min, self.step))
            self.x_vals = sorted(new_vals + self.x_vals)
        self.x_min = x_min
        self._calculate_y_vals()

    def set_x_max(self, x_max):
        if x_max < self.x_max:
            self.x_vals = [x for x in self.x_vals if x <= x_max]
        else:
            new_vals = list(range(self.x_max + 1, x_max + 1, self.step))
            self.x_vals.extend(new_vals)
        self.x_max = x_max
        self._calculate_y_vals()

    def set_y_min(self, y_min):
        self.y_min = y_min

    def set_y_max(self, y_max):
        self.y_max = y_max

    def set_params(self, params):
        self.params.update(params)
        self.br_calculator.set_params(self.params)
        self._calculate_y_vals()

    def add_channel(self, channel):
        if channel not in self.channels:
            self.channels.append(channel)
            self.y_vals[tuple(channel)] = []
            for x in self.x_vals:
                self.br_calculator.set_masses({"N1": x})
                self.y_vals[tuple(channel)].append(self.br_calculator.calculate(self.calculation_type, "N1", channel))

    def set_couplings(self, V1=None, V2=None, V3=None):
        if V1 is not None:
            self.params["V1"] = V1
        if V2 is not None:
            self.params["V2"] = V2
        if V3 is not None:
            self.params["V3"] = V3
        self.br_calculator.set_params(self.params)
        self._calculate_y_vals()

    def get_y_vals(self):
        return self.y_vals

# # Example usage:
# br_api = BR_API_Interface(1, 0, 10, 10, 1, list(range(1, 11)), "BR", [[14, -14, 14]])
# print(br_api.y_vals)
# br_api.add_channel([11, -11, 12])
# print(br_api.get_y_vals())
