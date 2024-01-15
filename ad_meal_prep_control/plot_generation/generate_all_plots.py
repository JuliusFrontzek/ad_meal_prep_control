import os

scenarios = ["1a", "1b", "2a", "2b", "2c"]

for scenario in scenarios:
    os.system(f"python3 plot_scenario_{scenario}.py 150 0")
