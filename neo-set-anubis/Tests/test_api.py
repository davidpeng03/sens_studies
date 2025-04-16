import unittest
from fastapi.testclient import TestClient
from App.api import app

class TestAPI(unittest.TestCase):

    def setUp(self):
        self.client = TestClient(app)

    def test_login_success(self):
        response = self.client.post("/token", data={"username": "johndoe", "password": "secret"})
        self.assertEqual(response.status_code, 200)
        self.assertIn("access_token", response.json())

    def test_login_failure(self):
        response = self.client.post("/token", data={"username": "johndoe", "password": "wrongpassword"})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {"detail": "Incorrect username or password"})

if __name__ == '__main__':
    unittest.main()
