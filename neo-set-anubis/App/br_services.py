# br_services.py
from App.utils.BR_API_interface import BR_API_Interface
from Core.Logger import CustomLogger
import math

br_api_instance = None

api_logger = CustomLogger('api', 'Logs/api.log').get_logger()

def initialize_br_api(x_min, y_min, x_max, y_max, step, x_vals, calculation_type, particles, model="HNL", params={"V1": 1, "V2": 1, "V3": 1}, masses={"N1": 1}):
    global br_api_instance
    api_logger.info(f"Initializing BR API with params: x_min={x_min}, y_min={y_min}, x_max={x_max}, y_max={y_max}, step={step}, x_vals={x_vals}, calculation_type={calculation_type}, particles={particles}, model={model}, params={params}, masses={masses}")
    br_api_instance = BR_API_Interface(x_min, y_min, x_max, y_max, step, x_vals, calculation_type, particles, model, params, masses)
    api_logger.info("BR API initialized successfully")

def set_x_min(x_min):
    if br_api_instance:
        api_logger.info(f"Setting x_min to {x_min}")
        br_api_instance.set_x_min(x_min)

def set_x_max(x_max):
    if br_api_instance:
        api_logger.info(f"Setting x_max to {x_max}")
        br_api_instance.set_x_max(x_max)

def set_y_min(y_min):
    if br_api_instance:
        api_logger.info(f"Setting y_min to {y_min}")
        br_api_instance.set_y_min(y_min)

def set_y_max(y_max):
    if br_api_instance:
        api_logger.info(f"Setting y_max to {y_max}")
        br_api_instance.set_y_max(y_max)

def set_params(params):
    if br_api_instance:
        api_logger.info(f"Setting params to {params}")
        br_api_instance.set_params(params)

def add_channel(channel):
    if br_api_instance:
        api_logger.info(f"Adding channel {channel}")
        br_api_instance.add_channel(channel)

def sanitize_y_vals(y_vals):
    sanitized = {}
    for key, values in y_vals.items():
        sanitized[str(key)] = [None if (isinstance(val, float) and math.isnan(val)) else val for val in values]
    return sanitized

def get_y_vals():
    if br_api_instance:
        y_vals = br_api_instance.get_y_vals()
        api_logger.info(f"Retrieved y_vals: {y_vals}")
        sanitized_y_vals = sanitize_y_vals(y_vals)
        api_logger.info(f"Sanitized y_vals: {sanitized_y_vals}")
        return sanitized_y_vals
    api_logger.warning("BR API instance is not initialized")
    return None


