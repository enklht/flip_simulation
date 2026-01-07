use crate::constants::{
    CELL_SIZE, DENSITY_STIFFNESS, DT, FLIP_RATIO, INCOMPLESSIBILITY_STEPS, OVER_RELAXATION,
};
use crate::get_acceleration;
use crate::grid::{CellType, MacGrid, b2f};
use crate::particle::Particle;

pub fn particles_to_grid(particles: &[Particle], grid: &mut MacGrid) {
    let iu = |i, j| i + j * (grid.nx + 1);
    let iv = |i, j| i + j * grid.nx;

    for p in particles {
        // cell type
        let i = (p.pos.x / CELL_SIZE).floor() as usize;
        let j = (p.pos.y / CELL_SIZE).floor() as usize;
        if i < grid.nx && j < grid.ny && grid.cell_type[i + j * grid.nx] == CellType::Air {
            grid.cell_type[i + j * grid.nx] = CellType::Water
        }

        // transfer u
        let x = p.pos.x / CELL_SIZE;
        let y = p.pos.y / CELL_SIZE - 0.5;

        let i0 = x.floor() as i32;
        let j0 = y.floor() as i32;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i > grid.nx as i32 || j >= grid.ny as i32 {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = iu(i as usize, j as usize);
                grid.u[idx] += p.vel.x * w;
                grid.wu[idx] += w;
            }
        }

        // transfer v
        let x = p.pos.x / CELL_SIZE - 0.5;
        let y = p.pos.y / CELL_SIZE;

        let i0 = x.floor() as i32;
        let j0 = y.floor() as i32;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as i32 || j > grid.ny as i32 {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = iv(i as usize, j as usize);
                grid.v[idx] += p.vel.y * w;
                grid.wv[idx] += w;
            }
        }
    }

    // normalize velocities
    for i in 0..grid.u.len() {
        if grid.wu[i] > 0. {
            grid.u[i] /= grid.wu[i];
        }
    }
    for i in 0..grid.v.len() {
        if grid.wv[i] > 0. {
            grid.v[i] /= grid.wv[i];
        }
    }
}

pub fn compute_particle_density(particles: &[Particle], grid: &mut MacGrid) {
    for p in particles {
        let x = p.pos.x / CELL_SIZE - 0.5;
        let y = p.pos.y / CELL_SIZE - 0.5;

        let i0 = x.floor() as i32;
        let j0 = y.floor() as i32;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as i32 || j >= grid.ny as i32 {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                grid.density[i as usize + j as usize * grid.nx] += w;
            }
        }
    }

    if grid.rest_density == 0. {
        let mut total_density = 0.;
        let mut water_cell_count = 0;

        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if grid.cell_type[i + j * grid.nx] == CellType::Water {
                    total_density += grid.density[i + j * grid.nx];
                    water_cell_count += 1;
                }
            }
        }

        if water_cell_count > 0 {
            grid.rest_density = total_density / water_cell_count as f32;
        }
    }
}

pub fn grid_to_particles(
    particles: &mut [Particle],
    grid: &MacGrid,
    u_prev: &[f32],
    v_prev: &[f32],
) {
    let iu = |i, j| i + j * (grid.nx + 1);
    let iv = |i, j| i + j * grid.nx;

    for p in particles {
        // transfer u
        let x = p.pos.x / CELL_SIZE;
        let y = p.pos.y / CELL_SIZE - 0.5;

        let i0 = x.floor() as i32;
        let j0 = y.floor() as i32;

        let mut ux = 0.;
        let mut ux_prev = 0.;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i > grid.nx as i32 || j >= grid.ny as i32 {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = iu(i as usize, j as usize);
                ux += grid.u[idx] * w;
                ux_prev += u_prev[idx] * w;
            }
        }

        // transfer v
        let x = p.pos.x / CELL_SIZE - 0.5;
        let y = p.pos.y / CELL_SIZE;

        let i0 = x.floor() as i32;
        let j0 = y.floor() as i32;

        let mut vy = 0.;
        let mut vy_prev = 0.;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as i32 || j > grid.ny as i32 {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = iv(i as usize, j as usize);
                vy += grid.v[idx] * w;
                vy_prev += v_prev[idx] * w;
            }
        }

        let pic = macroquad::prelude::vec2(ux, vy);
        let flip = p.vel + macroquad::prelude::vec2(ux - ux_prev, vy - vy_prev);
        p.vel = (1. - FLIP_RATIO) * pic + FLIP_RATIO * flip;
    }
}

pub fn apply_gravity(grid: &mut MacGrid) {
    let g = get_acceleration() * DT;

    grid.u.iter_mut().for_each(|u| *u += g.x);
    grid.v.iter_mut().for_each(|v| *v += g.y);
}

pub fn solve_incompressibility(grid: &mut MacGrid) {
    let nx = grid.nx;
    let ny = grid.ny;

    let is_wall = |i, j| grid.cell_type[i + j * nx] == CellType::Wall;
    let iu = |i, j| i + j * (nx + 1);
    let iv = |i, j| i + j * nx;

    for _ in 0..INCOMPLESSIBILITY_STEPS {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if grid.cell_type[i + j * nx] != CellType::Water {
                    continue;
                }

                let s_l = b2f(i > 0 && !is_wall(i - 1, j));
                let s_r = b2f(i < nx - 1 && !is_wall(i + 1, j));
                let s_t = b2f(j > 0 && !is_wall(i, j - 1));
                let s_b = b2f(j < ny - 1 && !is_wall(i, j + 1));

                let s = s_l + s_r + s_b + s_t;
                if s == 0. {
                    continue;
                }

                let mut d = OVER_RELAXATION
                    * (grid.u[iu(i + 1, j)] - grid.u[iu(i, j)] + grid.v[iv(i, j + 1)]
                        - grid.v[iv(i, j)]);

                d -= DENSITY_STIFFNESS * (grid.density[i + j * nx] - grid.rest_density);

                grid.u[iu(i, j)] += d * s_l / s;
                grid.u[iu(i + 1, j)] -= d * s_r / s;
                grid.v[iv(i, j)] += d * s_t / s;
                grid.v[iv(i, j + 1)] -= d * s_b / s;
            }
        }
    }
}
