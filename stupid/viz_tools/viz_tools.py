import numpy as np
import yt
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context
import os
import shutil
from multiprocessing import Pool
from functools import partial
import tqdm


class DataVisualizer:
    def __init__(self, file_prefix, num_cells):
        self.file_prefix = file_prefix
        self.num_cells = num_cells

    @staticmethod
    def load_binary_file(fname):
        data = np.fromfile(fname, dtype=np.double)

        a = data[0]
        particles = data[1:].reshape([-1, 3])

        return a, particles

    def load_snapshot(self, n):
        a, x = self.load_binary_file(self.file_prefix + '_x_{0}.bin'.format(n))
        _, p = self.load_binary_file(self.file_prefix + '_p_{0}.bin'.format(n))

        m = np.ones(x.shape[0])

        data = {
            'particle_position_x': x[:, 0],
            'particle_position_y': x[:, 1],
            'particle_position_z': x[:, 2],
            'particle_velocity_x': p[:, 0],
            'particle_velocity_y': p[:, 1],
            'particle_velocity_z': p[:, 2],
            'particle_mass': m
        }

        bbox = np.array([[0, self.num_cells], [0, self.num_cells], [0, self.num_cells]])

        ds = yt.load_particles(data, bbox=bbox)
        ds.current_time = a

        return ds

    def projection_plot(self, n, axis):
        self.ds = self.load_snapshot(n)

        p = yt.ProjectionPlot(self.ds, axis, ('deposit', 'all_density'))

        return p

    def animated_projection_plot(self, fname, snapshot_range, axis):
        p = self.projection_plot(snapshot_range[0], axis)
        fig = p[('deposit', 'all_density')].figure

        def animate(n):
            self.ds = self.load_snapshot(n)
            p._switch_ds(self.ds)

        animation = FuncAnimation(fig, animate, frames=snapshot_range)

        with rc_context({'mathtext.fontset': 'stix'}):
            animation.save(fname)

    def volume_plot(self, n, resolution=None):
        ds = self.load_snapshot(n)
        region = ds.r[::self.num_cells * 1j, ::self.num_cells * 1j, ::self.num_cells * 1j]

        density = region[('deposit', 'all_density')]
        data = {'density': density}

        bbox = np.array([[0, self.num_cells], [0, self.num_cells], [0, self.num_cells]])

        ds = yt.load_uniform_grid(data, density.shape, bbox=bbox)

        sc = yt.create_scene(ds, ('stream', 'density'))

        if resolution is not None:
            sc.camera.resolution = (resolution, resolution)

        return sc

    def _animated_volume_plot_worker(self, n, sigma_clip, resolution):
        sc = self.volume_plot(n, resolution)

        sc.save('tmp_renders/frame_{0}.png'.format(n), sigma_clip)

    def animated_volume_plot(self, fname, snapshot_range, sigma_clip=None, resolution=None, tqdm_notebook=False, parallel=False):
        os.makedirs('tmp_renders', exist_ok=True)

        if parallel:
            if tqdm_notebook:
                tqdm_type = tqdm.tqdm_notebook
            else:
                tqdm_type = tqdm.tqdm

            with Pool() as pool:
                list(tqdm_type(pool.imap(
                    partial(self._animated_volume_plot_worker, sigma_clip=sigma_clip, resolution=resolution),
                    snapshot_range), total=len(snapshot_range)))
        else:
            for n in snapshot_range:
                sc = self.volume_plot(n, resolution)

            sc.save('tmp_renders/frame_{0}.png'.format(n), sigma_clip)

        os.system('ffmpeg -f image2 -r 15 -i tmp_renders/frame_%d.png -vcodec h264 -y {0}'.format(fname))

        shutil.rmtree('tmp_renders')