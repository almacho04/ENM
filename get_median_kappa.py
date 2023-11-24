import json
# file_path = "./DB_1/sENM10.json"
# file_path = "./DB_1/sdENM.json"
# file_path = "./DB_1/dENM.json"
file_path = "./DB_1/sENM13.json"
with open(file_path, "r") as f:
    data = json.load(f)

kappa_values = [item["kappa"] for item in data]
kappa_values.sort()

if len(kappa_values) % 2 == 0:
    median = (kappa_values[len(kappa_values)//2] + kappa_values[len(kappa_values)//2 - 1]) / 2
else:
    median = kappa_values[len(kappa_values)//2]

print(f"Median kappa value: {median}")