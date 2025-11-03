import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
N = 1000
timesteps = 50
beta = 0.02           # transmission probability per contact
clear_prob = 0.05     # clearance probability per timestep
vacc_coverage = 0.7
vaccine_efficacy = 0.9

# Initialize population
vaccinated = np.random.rand(N) < vacc_coverage
infected = np.zeros(N, dtype=bool)
infected[np.random.randint(0, N//20)] = True  # small seed

prevalence = []

for t in range(timesteps):
    for i in range(N):
        if infected[i]:
            partners = np.random.choice(N, size=5, replace=False)
            for p in partners:
                if not infected[p]:
                    p_eff = beta * (1 - vaccine_efficacy if vaccinated[p] else 1)
                    if np.random.rand() < p_eff:
                        infected[p] = True
            # Clearance
            if np.random.rand() < clear_prob:
                infected[i] = False
    prevalence.append(infected.mean())

plt.plot(prevalence)
plt.xlabel('Time step')
plt.ylabel('Prevalence')
plt.title('Preliminary HPV-ABM Simulation')
plt.show()
