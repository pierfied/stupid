import numpy as np
import yt
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context


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

        data = {
            'particle_position_x': x[:, 0],
            'particle_position_y': x[:, 1],
            'particle_position_z': x[:, 2],
            'particle_velocity_x': p[:, 0],
            'particle_velocity_y': p[:, 1],
            'particle_velocity_z': p[:, 2]
        }

        bbox = np.array([[0, self.num_cells], [0, self.num_cells], [0, self.num_cells]])

        ds = yt.load_particles(data, bbox=bbox)
        ds.current_time = a

        return ds

    def projection_plot(self, n, axis):
        self.ds = self.load_snapshot(n)

        p = yt.ProjectionPlot(self.ds, axis, ('deposit', 'all_count'))

        return p

    def animated_projection_plot(self, fname, snapshot_range, axis):
        p = self.projection_plot(snapshot_range[0], axis)
        fig = p[('deposit', 'all_count')].figure

        def animate(n):
            self.ds = self.load_snapshot(n)
            p._switch_ds(self.ds)

        animation = FuncAnimation(fig, animate, frames=snapshot_range)

        with rc_context({'mathtext.fontset': 'stix'}):
            animation.save(fname)
