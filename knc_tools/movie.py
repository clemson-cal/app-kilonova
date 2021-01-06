import argparse
import matplotlib.pyplot as plt
import numpy as np
from knc_tools import products


light_speed = 3e10


def pcolormesh_blocks(ax, pcm, data_range):
    vmin = data_range[0] if data_range[0] is not None else min(C.min() for X, Z, C in pcm)
    vmax = data_range[1] if data_range[1] is not None else max(C.max() for X, Z, C in pcm)

    for X, Z, C in pcm:
        C[C != C] = vmin
        cm = ax.pcolormesh(X, Z, C, vmin=vmin, vmax=vmax, edgecolors='none')
    return cm


def plot(fig, filename):
    spec = fig.add_gridspec(2, 3, height_ratios=[19, 1])
    ax1  = fig.add_subplot(spec[0, 0])
    ax2  = fig.add_subplot(spec[0, 1])
    ax3  = fig.add_subplot(spec[0, 2])
    cax1 = fig.add_subplot(spec[1, 0])
    cax2 = fig.add_subplot(spec[1, 1])
    cax3 = fig.add_subplot(spec[1, 2])

    prods = products.Products(filename)
    tstart = prods.config['control']['start_time']
    rstart = 2e8
    router = rstart + (prods.time - tstart) * light_speed
    rho_in = 1e13 * (router / rstart)**-2
    pre_in = rho_in * 1e-3

    for ax in [ax1, ax2, ax3]:
        ax.set_aspect('equal')
        ax.set_xlim(0.0, router)
        ax.set_ylim(0.0, router)
        ax.set_xlabel(r'$x [\rm{cm}]$')

    pcm_d = prods.pcolormesh_data('rho', log=True)
    pcm_p = prods.pcolormesh_data('pre', log=True)
    pcm_u = prods.pcolormesh_data('ur',  log=False)

    cm1 = pcolormesh_blocks(ax1, pcm_d, [np.log10(rho_in), np.log10(rho_in * 1e5)])
    cm2 = pcolormesh_blocks(ax2, pcm_p, [np.log10(pre_in), np.log10(pre_in * 1e5)])
    cm3 = pcolormesh_blocks(ax3, pcm_u, [0.0, 10.0])

    ax1.set_title('Density')
    ax2.set_title('Pressure')
    ax3.set_title(r'Radial $\Gamma \beta$')
    fig.colorbar(cm1, cax1, orientation='horizontal')
    fig.colorbar(cm2, cax2, orientation='horizontal')
    fig.colorbar(cm3, cax3, orientation='horizontal')
    fig.suptitle(r'$t = {:.02f} s$'.format(prods.time))


def make_movie(fig, filenames, output):
    from matplotlib.animation import FFMpegWriter
    writer = FFMpegWriter(fps=15)
    with writer.saving(fig, output, dpi=200):
        for filename in filenames:
            print(filename)
            plot(fig, filename)
            writer.grab_frame()
            fig.clf()
    print('writing', output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', default='output.mp4')
    args = parser.parse_args()

    fig = plt.figure(figsize=[16, 7])
    fig.subplots_adjust(top=0.998, right=0.95, left=0.05, bottom=0.05)

    if len(args.filenames) > 1:
        make_movie(fig, sorted(args.filenames), args.output)
    else:
        plot(fig, args.filenames[0])
        plt.show()
