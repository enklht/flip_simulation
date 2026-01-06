use macroquad::{miniquad::window::screen_size, prelude::*};
use std::collections::HashMap;

const DT: f32 = 1. / 24.;
const GRAVITY: Vec2 = vec2(0., 300.);
const RADIUS: f32 = 5.;
const CELL_SIZE: f32 = 10.;

const OVER_RELAXATION: f32 = 1.9;
const DENSITY_STIFFNESS: f32 = 1.;
const FLIP_RATIO: f32 = 0.9;
const SEPARATION_STEPS: usize = 3;
const INCOMPLESSIBILITY_STEPS: usize = 5;

struct Particle {
    pos: Vec2,
    vel: Vec2,
}

#[derive(Debug, Clone, PartialEq)]
enum CellType {
    Air,
    Water,
    Wall,
}

struct MacGrid {
    nx: usize,
    ny: usize,
    u: Vec<f32>,
    v: Vec<f32>,
    wu: Vec<f32>,
    wv: Vec<f32>,
    cell_type: Vec<CellType>,
    density: Vec<f32>,
    rest_density: f32,
}

impl MacGrid {
    fn new(w: f32, h: f32) -> Self {
        let nx = (w / CELL_SIZE).floor() as usize;
        let ny = (h / CELL_SIZE).floor() as usize;
        let mut cell_type = vec![CellType::Air; nx * ny];

        for j in 0..ny {
            cell_type[j * nx] = CellType::Wall;
            cell_type[j * nx + nx - 1] = CellType::Wall;
        }
        for i in 0..nx {
            cell_type[i] = CellType::Wall;
            cell_type[(ny - 1) * nx + i] = CellType::Wall;
        }

        Self {
            nx,
            ny,
            u: vec![0.; (nx + 1) * ny],
            v: vec![0.; nx * (ny + 1)],
            wu: vec![0.; (nx + 1) * ny],
            wv: vec![0.; nx * (ny + 1)],
            cell_type,
            density: vec![0.; nx * ny],
            rest_density: 0.,
        }
    }

    fn clear(&mut self) {
        self.u.fill(0.);
        self.v.fill(0.);
        self.wu.fill(0.);
        self.wv.fill(0.);
        self.density.fill(0.);
        self.cell_type.fill(CellType::Air);

        let nx = self.nx;
        let ny = self.ny;

        for j in 0..ny {
            self.cell_type[j * nx] = CellType::Wall;
            self.cell_type[j * nx + nx - 1] = CellType::Wall;
        }
        for i in 0..nx {
            self.cell_type[i] = CellType::Wall;
            self.cell_type[(ny - 1) * nx + i] = CellType::Wall;
        }
    }
}

fn integrate_particles(particles: &mut [Particle]) {
    for p in particles {
        p.vel += GRAVITY * DT;
        p.pos += p.vel * DT;
    }
}

fn push_particles_apart(particles: &mut [Particle]) {
    let mut map = HashMap::<(usize, usize), Vec<usize>>::new();

    for (i, p) in particles.iter().enumerate() {
        let cell = (
            (p.pos.x / CELL_SIZE).floor() as usize,
            (p.pos.y / CELL_SIZE).floor() as usize,
        );
        map.entry(cell).or_default().push(i);
    }

    for _ in 0..SEPARATION_STEPS {
        for i in 0..particles.len() {
            let cell = (
                (particles[i].pos.x / CELL_SIZE).floor() as usize,
                (particles[i].pos.y / CELL_SIZE).floor() as usize,
            );

            for x in cell.0.saturating_sub(1)..=cell.0 + 1 {
                for y in cell.1.saturating_sub(1)..=cell.1 + 1 {
                    let Some(neighbours) = map.get(&(x, y)) else {
                        continue;
                    };

                    for &j in neighbours {
                        if j <= i {
                            continue;
                        }

                        let delta = particles[j].pos - particles[i].pos;
                        let d = delta.length();
                        if d == 0. || d >= 2. * RADIUS {
                            continue;
                        }

                        let s = 0.5 * (2. * RADIUS - d) / d;

                        particles[i].pos -= delta * s;
                        particles[j].pos += delta * s;
                    }
                }
            }
        }
    }
}

fn handle_particle_collisions(particles: &mut [Particle], grid: &MacGrid) {
    for p in particles {
        let min_water_x = CELL_SIZE + RADIUS;
        let max_water_x = (grid.nx - 1) as f32 * CELL_SIZE - RADIUS;
        let min_water_y = CELL_SIZE + RADIUS;
        let max_water_y = (grid.ny - 1) as f32 * CELL_SIZE - RADIUS;

        if p.pos.x < min_water_x {
            p.pos.x = min_water_x;
            p.vel.x = 0.;
        }
        if p.pos.x > max_water_x {
            p.pos.x = max_water_x;
            p.vel.x = 0.;
        }
        if p.pos.y < min_water_y {
            p.pos.y = min_water_y;
            p.vel.y = 0.;
        }
        if p.pos.y > max_water_y {
            p.pos.y = max_water_y;
            p.vel.y = 0.;
        }
    }
}

fn particles_to_grid(particles: &[Particle], grid: &mut MacGrid) {
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

fn compute_particle_density(particles: &[Particle], grid: &mut MacGrid) {
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

fn grid_to_particles(particles: &mut [Particle], grid: &MacGrid, u_prev: &[f32], v_prev: &[f32]) {
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

        let pic = vec2(ux, vy);
        let delta = vec2(ux - ux_prev, vy - vy_prev);
        p.vel = (1. - FLIP_RATIO) * pic + FLIP_RATIO * (p.vel + delta);
    }
}

fn spawn_particles(n: usize) -> Vec<Particle> {
    let mut particles = Vec::with_capacity(n);

    for i in 0..n {
        let x = 100. + (i as f32 % 50.) * CELL_SIZE;
        let y = 50. + (i as f32 / 50.).floor() * CELL_SIZE;

        particles.push(Particle {
            pos: vec2(x, y),
            vel: vec2(0., 0.),
        })
    }

    particles
}

#[inline]
fn b2f(b: bool) -> f32 {
    if b { 1. } else { 0. }
}

fn solve_incompressibility(grid: &mut MacGrid) {
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

fn render(particles: &[Particle], grid: &MacGrid) {
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            if grid.cell_type[i + j * grid.nx] == CellType::Wall {
                let x = i as f32 * CELL_SIZE;
                let y = j as f32 * CELL_SIZE;
                let color = Color::new(0.5, 0.5, 0.5, 1.0);
                draw_rectangle(x, y, CELL_SIZE, CELL_SIZE, color);
            }
        }
    }

    for p in particles {
        draw_circle(p.pos.x, p.pos.y, RADIUS, SKYBLUE);
    }
}

#[macroquad::main("FLIP")]
async fn main() {
    let (mut w, mut h) = screen_size();

    let mut particles = spawn_particles(10000);
    let mut grid = MacGrid::new(w, h);

    loop {
        clear_background(BLACK);

        let new_size = screen_size();
        if (w, h) != new_size {
            (w, h) = new_size;
            grid = MacGrid::new(w, h);
        } else {
            grid.clear();
        }
        grid.clear();

        integrate_particles(&mut particles);
        push_particles_apart(&mut particles);
        handle_particle_collisions(&mut particles, &grid);

        particles_to_grid(&particles, &mut grid);

        compute_particle_density(&particles, &mut grid);

        let u_prev = grid.u.clone();
        let v_prev = grid.v.clone();

        solve_incompressibility(&mut grid);

        grid_to_particles(&mut particles, &grid, &u_prev, &v_prev);

        render(&particles, &grid);

        next_frame().await
    }
}
