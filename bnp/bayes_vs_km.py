import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from itertools import product
import seaborn as sns
import matplotlib.colors as colors

# generate censored data given two scipy distributions
def generate_data(n, f_dist, g_dist):
    f_sample = f_dist.rvs(n)
    g_sample = g_dist.rvs(n)

    observations = np.minimum(f_sample, g_sample)

    deltas = f_sample < g_sample

    return observations, deltas.astype(int)

# counting process: number of non-censored observations
def countingN(time, obs, deltas):
    valid_obs = (obs <= time)
    valid_obs.astype(int)

    prod = valid_obs * deltas

    return np.sum(prod)

# counting process: individuals at risk
def countingY(time, obs, deltas):
    valid_obs = (obs >= time)
    return np.sum(valid_obs)

# compute KM estimator on a grid
def KaplanMeierEstimator(time_grid, obs, deltas):
    prod = 1
    KME = [prod]
    for i in range(0, len(time_grid)-1):
        t1 = time_grid[i]
        t2 = time_grid[i+1]

        delta_n = countingN(t2, obs, deltas) - countingN(t1, obs, deltas)
        Y = countingY(t1, obs, deltas)
        if Y > 0:
            prod = prod * (1 - delta_n / Y)
        else:
            prod = 0

        KME.append(prod)
    
    return np.array(KME)

# compute DE estimator on a grid
def DirichletEstimator(time_grid, obs, deltas, prior, c = 1):
    prod = 1
    KME = [prod]
    for i in range(0, len(time_grid)-1):
        t1 = time_grid[i]
        t2 = time_grid[i+1]

        alpha1 = (prior.cdf(t2) - prior.cdf(t1))
        alpha2 = 1 - prior.cdf(t1)

        if alpha2 < alpha1:
            print('ERROR 1')

        delta_n = countingN(t2, obs, deltas) - countingN(t1, obs, deltas)
        Y = countingY(t1, obs, deltas)
        if Y < delta_n:
            print(delta_n, Y)
            print('ERROR 2')

        prod = prod * (1 - (c*alpha1 + delta_n) / (c*alpha2 + Y))
            
        KME.append(prod)
    
    return np.array(KME)

# calculate ISE using trapezium rule
def ISE(S, S_est, time_grid, lamb):
    out = 0
    for i in range(0, len(time_grid)-1):
        t2 = time_grid[i+1]
        t1 = time_grid[i]

        f_left = (S[i] - S_est[i])**2
        f_right = (S[i+1] - S_est[i+1])**2

        out += (f_left + f_right)/2 * (t2 - t1)
    return out

# concatenate the grid and the observed times
def full_grid(time_grid, x_obs):
    full_grid = np.concatenate((time_grid, x_obs))
    return np.sort(full_grid)


# ------------------------------------------------------------------------

# This part of the code generates the matrices and saves them


lambdas = [.5, 1, 2, 3, 5, 10, 15]
alphas = [1, 5, 10, 25, 50, 100]

n_avg = 100

lam_len = len(lambdas)
alph_len = len(alphas)

time_grid0 = np.linspace(0, 2.5, 500)
dt = time_grid0[1] - time_grid0[0]

f = scipy.stats.expon(scale = 1 / 2)
g = scipy.stats.expon(scale = 1 / 2)

indexes = list(product(range(lam_len), range(alph_len)))

for n_obs in [10, 25, 75, 150]:

    results = np.zeros((lam_len, alph_len))

    for i, j in indexes:
        lamb = lambdas[i]
        alpha = alphas[j]

        prior = scipy.stats.expon(scale = 1 / lamb)

        in_loop_results = []

        for k in range(0, n_avg):
            obs, deltas = generate_data(n_obs, f, g)
            time_grid = full_grid(time_grid0, obs)
            survivor_de = DirichletEstimator(time_grid, obs, deltas, prior, c = alpha)
            survior_real = f.sf(time_grid)
            in_loop_results.append(ISE(survior_real, survivor_de, time_grid, 2))

        results[i, j] = np.mean(in_loop_results)
        print('Done with index: i = ', i, ' ; j = ', j, ' for n = ', n_obs)

    print(results)

    in_loop_results = []
    for k in range(0, n_avg):
        obs, deltas = generate_data(n_obs, f, g)
        time_grid = full_grid(time_grid0, obs)
        survivor_km = KaplanMeierEstimator(time_grid, obs, deltas)
        in_loop_results.append(ISE(survior_real, survivor_km, time_grid, 2))

    ke_val = np.mean(in_loop_results)

    mat = ke_val - results

    file_name = 'matrix_n_' + str(n_obs)
    np.save(file_name, mat)

# -------------------------------------

# This part of the code generates the figure and saves it

n_list = [10, 25, 75, 150]

mats = []
max_val = 0
min_val = 0

for n in n_list:
    mat = np.load('matrix_n_' + str(n) + '.npy')
    mats.append(mat)
    max_val = max(max_val, mat.max())
    min_val = min(min_val, mat.min())
    

# Two slope norm for color bar
colors_positive = plt.cm.Reds(np.linspace(1, .1, 256))
colors_negative = plt.cm.Blues(np.linspace(.1, 1, 256))
all_colors = np.vstack((colors_positive, colors_negative))
color_map = colors.LinearSegmentedColormap.from_list('color_map', all_colors)

divnorm = colors.TwoSlopeNorm(vcenter = 0, vmin=min_val, vmax=max_val)

# initalise subplots
fig, ax = plt.subplots(nrows=2, ncols = 2)

fig.set_figheight(12)
fig.set_figwidth(14)

#increase space between plots
plt.subplots_adjust(hspace = 0.3)
fig.suptitle('Comparison between Bayesian and KM Estimators')

# axes
axes = [[0, 0], [0, 1], [1, 0], [1, 1]]

mat_num = 0
for i, j in axes:
    n = n_list[mat_num]
    mat = mats[mat_num]
    mat_num += 1
    sns.heatmap(mat, ax = ax[i, j], robust = True, xticklabels = alphas, yticklabels = lambdas, cmap = color_map, cbar = False, norm = divnorm)
    ax[i, j].set_xlabel(r'c')
    ax[i, j].set_ylabel(r'$\lambda_{prior}$')
    ax[i, j].tick_params(left = False, bottom = False)
    # set lines between 
    ax[i, j].hlines(range(1, 7), *ax[i, j].get_xlim(), linewidth = 3)
    ax[i, j].vlines(range(1, 7), *ax[i, j].get_ylim(), linewidth = 3)
    ax[i, j].title.set_text('n = ' + str(n))

fig.colorbar(ax[0, 0].get_children()[0], ax = ax.ravel().tolist(), orientation = 'horizontal', label = r'AVG($ISE_{KM} - ISE_{DP}$)', fraction = 0.04, pad = 0.1)

plt.savefig('final_graph', dpi = 800)