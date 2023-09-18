
import numpy as np

# Data source for Rn-224 chain:
# https://www.epa.gov/sites/default/files/2018-12/documents/radon-22-decay-chain.pdf

second = 1.0
minute = 60  * second
hour   = 60  * minute
day    = 24  * hour
year   = 365 * day

# Type of decay (alpha, or beta)
rad = ["a", "a", "b", "b", "a", "b", "b", "a"]

# Decay time constants for radon 224 chain
# Rn-224 -> Po-218 -> Pb-214 -> Bi-214 -> Po-214 -> Pb-210 -> Bi-210 -> Po-210 >- Pb-206
tau = np.array([3.8*day, 3*minute, 28.6*minute, 19.7*minute, 0.0002*second, 22*year, 5*day, 138*day])

# Isotope populations
pop = [0 for _ in tau]
pop.append(0) # add an extra population on the end for Pb-206

# Room with only Radon
pop[0] = 1

# Calculate dN for a 1 second time period
decayRate = np.power(0.5, second/tau)

# Create store arrays for alpha and beta radiation
alphaRad = []
betaRad = []


# Simulate decay chain
simtime = 0
calculating = True
while simtime < (tau[0]*2 + 60 * minute):

    # Run decay process from radon-224 (allow extra time for sim to stabilise)
    # For your case, Matthias, the balloon is currently collecting decay products.

    newpop = [s*t for (s,t) in zip(decayRate, pop)] # non decayed fraction for each species
    dN = [p-np for (p,np) in zip(pop, newpop)] # change for each species
    alpha = sum([n for (n,rad) in zip(dN,rad) if rad == "a" ]) # alpha radiation
    beta = sum([n for (n,rad) in zip(dN,rad) if rad == "b" ]) # beta radiation
    dN.insert(0,0)
    pop = [n+d for (n,d) in zip(newpop, dN)] # calculate new populations
    simtime += 1


    # Let Rn-222 form decay products until simulation is 'stable'
    if simtime > tau[0]*2:
        # From here, remove the balloon from collecting decay prods.
        # This is the same as removing all Rn-224 from the room.
        pop[0] = 0
        # Store radiation measurements
        alphaRad.append(alpha)
        betaRad.append(beta)

# Trim the first element to ignore edge effects from instantanous randon removal
del alphaRad[0]
del betaRad[0]


# Plot the results
import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)

maxRad = max(max(alphaRad), max(betaRad))

t = [x/minute for x, _ in enumerate(alphaRad)]


# show radiation on linear plot
fig, ax = plt.subplots()
betaline, = ax.plot(t, betaRad/maxRad, label="Beta Radiation")
alphaline, = ax.plot(t, alphaRad/maxRad, label="Alpha Radiation")
plt.xlabel("Time [mins]")
plt.ylabel("Radiation")
ax.legend(handles=[betaline, alphaline])


logbeta = np.log2(betaRad/maxRad)
logalpha = np.log2(alphaRad/maxRad)

# show radiation on log scale
fig, ax = plt.subplots()
betaline,  = ax.plot(t, logbeta, label="Beta Radiation")
alphaline, = ax.plot(t, logalpha, label="Alpha Radiation")
plt.xlabel("Time [mins]")
plt.ylabel("Log2 Radiation")
ax.legend(handles=[betaline, alphaline])

# estimate the instantanous time constant (I think this is correct...)
fig, ax = plt.subplots()
betaline,  = ax.plot(t, -1/(np.gradient(logbeta)*minute), label="Beta Radiation")
alphaline, = ax.plot(t, -1/(np.gradient(logalpha)*minute), label="Alpha Radiation")
plt.xlabel("Time [mins]")
plt.ylabel("Instantaneous half life [mins]")
ax.legend(handles=[betaline, alphaline])
ax.set_ylim(0, 80)

plt.show()
