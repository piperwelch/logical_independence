from utils import * 
import matplotlib.pyplot as plt 
import glob 
import numpy as np 

eval_method = "simult"
seed = 0 
id = 19786
output = 23
fig, ax = plt.subplots(4, 4, sharex=True, sharey=True, figsize = (10, 6))
ax = ax.flatten()
files = glob.glob(f"dump.{eval_method}_eval_seed{seed}_id{id}*")
files.sort()
for index, file in enumerate(files):
    index+=1
    condition = file.split("_")[-1]
    data_dict = read_particle_data(file)

    # print(data_dict[output])
    xs, ys = [], []
    for pt in data_dict[output]:
        xs.append(pt[0])
        ys.append(pt[1])
    xs = np.array(xs) - np.mean(xs)
    # quit()
    # ax[index].plot(xs)
    # ax[index].set_title(condition)
    fft_xs = np.fft.rfft(xs)[28:50]
    fft_freqs = np.abs(np.fft.rfftfreq(len(xs), d=0.005))[28:50]
    print(fft_freqs)
    # Plot FFT
    ax[index].plot(fft_freqs, np.abs(fft_xs))
    ax[index].set_title(f'{condition} - FFT')
# plt.axis('off')

plt.tight_layout()
plt.savefig("plot")
quit()