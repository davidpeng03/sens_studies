from pydantic import BaseModel

class UserInDB(BaseModel):
    username: str
    full_name: str
    email: str
    hashed_password: str
    disabled: bool = False
