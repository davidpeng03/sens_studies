from fastapi import FastAPI, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from .schemas import Token, User
from .services import fake_users_db, fake_hash_password
from .auth import get_current_active_user
from .models import UserInDB
from . import br_services
from pydantic import BaseModel
from typing import List, Dict, Any
from Core.Logger import CustomLogger

app = FastAPI()

api_logger = CustomLogger('api', 'Logs/api.log').get_logger()

@app.post("/token", response_model=Token)
async def login(form_data: OAuth2PasswordRequestForm = Depends()):
    user_dict = fake_users_db.get(form_data.username)
    if not user_dict:
        raise HTTPException(status_code=400, detail="Incorrect username or password")
    user = UserInDB(**user_dict)
    hashed_password = fake_hash_password(form_data.password)
    if not hashed_password == user.hashed_password:
        raise HTTPException(status_code=400, detail="Incorrect username or password")

    return {"access_token": user.username, "token_type": "bearer"}

@app.get("/users/me", response_model=User)
async def read_users_me(current_user: User = Depends(get_current_active_user)):
    return current_user

@app.get("/simulate")
async def simulate(step: int, current_user: User = Depends(get_current_active_user)):
    return {"message": f"Running simulation step {step} for user {current_user.username}"}

class BRConfig(BaseModel):
    x_min: float
    y_min: float
    x_max: float
    y_max: float
    step: float
    x_vals: List[float]
    calculation_type: str
    particles: List[List[int]]
    model: str = "HNL"
    params: Dict[str, Any] = {"V1": 1, "V2": 1, "V3": 1}
    masses: Dict[str, int] = {"N1": 1}
    
@app.post("/br/init")
async def initialize_br(config: BRConfig, current_user: User = Depends(get_current_active_user)):
    br_services.initialize_br_api(
        config.x_min, config.y_min, config.x_max, config.y_max, config.step, config.x_vals,
        config.calculation_type, config.particles, config.model, config.params, config.masses
    )
    return {"message": "BR API initialized successfully"}

@app.post("/br/set_x_min")
async def set_x_min(x_min: int, current_user: User = Depends(get_current_active_user)):
    br_services.set_x_min(x_min)
    return {"message": "x_min set successfully"}

@app.post("/br/set_x_max")
async def set_x_max(x_max: int, current_user: User = Depends(get_current_active_user)):
    br_services.set_x_max(x_max)
    return {"message": "x_max set successfully"}

@app.post("/br/set_y_min")
async def set_y_min(y_min: int, current_user: User = Depends(get_current_active_user)):
    br_services.set_y_min(y_min)
    return {"message": "y_min set successfully"}

@app.post("/br/set_y_max")
async def set_y_max(y_max: int, current_user: User = Depends(get_current_active_user)):
    br_services.set_y_max(y_max)
    return {"message": "y_max set successfully"}

@app.post("/br/set_params")
async def set_params(params: Dict[str, Any], current_user: User = Depends(get_current_active_user)):
    br_services.set_params(params)
    return {"message": "Parameters set successfully"}

@app.post("/br/add_channel")
async def add_channel(channel: List[int], current_user: User = Depends(get_current_active_user)):
    br_services.add_channel(channel)
    return {"message": "Channel added successfully"}

@app.get("/br/get_y_vals")
async def get_y_vals(current_user: User = Depends(get_current_active_user)):
    try:
        y_vals = br_services.get_y_vals()
        api_logger.info(f"Returning y_vals: {y_vals}")
        return {"y_vals": y_vals}
    except Exception as e:
        api_logger.error(f"Get Y Vals error: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")
