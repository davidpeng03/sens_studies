from App.models import UserInDB

fake_users_db = {
    "theo": {
        "username": "theo",
        "full_name": "Theo Reymermier",
        "email": "theo.reymermier@cern.ch",
        "hashed_password": "fakehashedanubis",
        "disabled": False,
    },
    "johndoe": {
        "username": "johndoe",
        "full_name": "John Doe",
        "email": "johndoe@example.com",
        "hashed_password": "fakehashedsecret",
        "disabled": False,
    }
}

def fake_hash_password(password: str):
    return "fakehashed" + password

def get_user(db, username: str):
    if username in db:
        user_dict = db[username]
        return UserInDB(**user_dict)

def fake_decode_token(token):
    user = get_user(fake_users_db, token)
    return user

