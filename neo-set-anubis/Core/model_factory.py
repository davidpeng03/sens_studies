from Core.Ifactory import IModelFactory
from Core.HNL.hnl_factory import HNLFactory
from Core.DarkPhoton.darkphoton_factory import DarkPhotonFactory
from Core.GeneralBSM.general_factory import GeneralBSMFactory
class ModelFactory:
    @staticmethod        
    def get_factory(model_name : str) -> IModelFactory:
        if model_name == "HNL":
            return HNLFactory()
        elif model_name == "DarkPhoton":
            return DarkPhotonFactory()
        elif model_name == "GeneralBSM":
            return GeneralBSMFactory()
        else:
            raise ValueError(f"Unknown model type: {model_name}")