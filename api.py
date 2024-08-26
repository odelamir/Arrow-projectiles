from flask import Flask, request, jsonify
from flask_cors import CORS
import multiprocessing
import runpy

def run_another():
    runpy.run_path('run_c.py')

app = Flask(__name__)
cors = CORS(app)

users = [
    {"email": "fifa@gmail.com", "password": "1234"},
    {"email": "Elioo@gmail.com", "password": "8989"},
    {"email": "Amit@gmail.com", "password": "4545"},
    {"email": "odel@gmail.com", "password": "1234"},
    {"email": "Or@gmail.com", "password": "5566"}
]

@app.route("/map", methods=["GET"])
def Mark_threats():
  remaining_lines=0
  local_yrut_x = 0
  local_yrut_y =0
  local_yrut_z = 0
  current_time=0
  try:  
    with open('data.txt', 'r') as file:
        lines = file.readlines()
        if len(lines) >= 4:
            local_yrut_x = float(lines[0].strip())
            local_yrut_y = float(lines[1].strip())
            local_yrut_z = float(lines[2].strip())
            current_time=lines[3].strip()
            remaining_lines = lines[4:]
            if(remaining_lines):
             with open('data.txt', 'w') as file:
              file.writelines(remaining_lines)
            else:
                open('data.txt', 'w').close()             
  except FileNotFoundError:
   pass
  return jsonify({"points_now": {"x": local_yrut_x, "y": local_yrut_y, "z": local_yrut_z, "current_time": current_time}})


@app.route("/login", methods=["POST"])
def login():
    data = request.get_json()
    email = data.get("myemail")
    password = data.get("mypassword")
    print(email, password)
    for user in users:
        if user["email"] == email and user["password"] == password:
            return {"user": user}
    return "not found"


@app.route("/signup", methods=["POST"])
def signup():
    data = request.get_json()
    id = data.get("id")
    name = data.get("name")
    email = data.get("email")
    password = data.get("password")

    if not id or not name or not email or not password:
        return jsonify({"message": "Invalid input"}), 400

    if any(user["email"] == email for user in users):
        return jsonify({"message": "User already exists"}), 409

    new_user = {"id": id, "name": name, "email": email, "password": password}
    users.append(new_user)
    return jsonify({"message": "Add Successfully"}), 201

if __name__ == "__main__":
    multiprocessing.Process(target=run_another).start()
    app.run()