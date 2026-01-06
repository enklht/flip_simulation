use macroquad::prelude::*;

const DT: f32 = 1. / 24.;
const GRAVITY: Vec2 = vec2(0., 500.);
const RADIUS: f32 = 5.;
const CELL_SIZE: f32 = 12.;

const OVER_RELAXATION: f32 = 1.9;
const DENSITY_STIFFNESS: f32 = 3.;
const FLIP_RATIO: f32 = 0.98;

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

    fn iu(&self, i: usize, j: usize) -> usize {
        i + j * (self.nx + 1)
    }

    fn iv(&self, i: usize, j: usize) -> usize {
        i + j * self.nx
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

fn particles_to_grid(particles: &[Particle], grid: &mut MacGrid) {
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

        let i0 = x.floor() as isize;
        let j0 = y.floor() as isize;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i > grid.nx as isize || j >= grid.ny as isize {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = grid.iu(i as usize, j as usize);
                grid.u[idx] += p.vel.x * w;
                grid.wu[idx] += w;
            }
        }

        // transfer v
        let x = p.pos.x / CELL_SIZE - 0.5;
        let y = p.pos.y / CELL_SIZE;

        let i0 = x.floor() as isize;
        let j0 = y.floor() as isize;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as isize || j > grid.ny as isize {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                let idx = grid.iv(i as usize, j as usize);
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

        let i0 = x.floor() as isize;
        let j0 = y.floor() as isize;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as isize || j >= grid.ny as isize {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                grid.density[i as usize + j as usize * grid.nx] += w;
            }
        }
    }
}

fn grid_to_particles(particles: &mut [Particle], grid: &MacGrid, u_prev: &[f32], v_prev: &[f32]) {
    for p in particles {
        // transfer u
        let x = p.pos.x / CELL_SIZE;
        let y = p.pos.y / CELL_SIZE - 0.5;

        let i0 = x.floor() as isize;
        let j0 = y.floor() as isize;

        let mut ux = 0.;
        let mut ux_prev = 0.;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i > grid.nx as isize || j >= grid.ny as isize {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                ux += grid.u[grid.iu(i as usize, j as usize)] * w;
                ux_prev += u_prev[grid.iu(i as usize, j as usize)] * w;
            }
        }

        // transfer v
        let x = p.pos.x / CELL_SIZE - 0.5;
        let y = p.pos.y / CELL_SIZE;

        let i0 = x.floor() as isize;
        let j0 = y.floor() as isize;

        let mut vy = 0.;
        let mut vy_prev = 0.;

        for di in 0..=1 {
            for dj in 0..=1 {
                let i = i0 + di;
                let j = j0 + dj;

                if i < 0 || j < 0 || i >= grid.nx as isize || j > grid.ny as isize {
                    continue;
                }

                let wx = 1. - (x - i as f32).abs();
                let wy = 1. - (y - j as f32).abs();
                let w = wx * wy;

                vy += grid.v[grid.iv(i as usize, j as usize)] * w;
                vy_prev += v_prev[grid.iv(i as usize, j as usize)] * w;
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
        let x = 100. + (i as f32 % 50.) * 12.;
        let y = 50. + (i as f32 / 50.).floor() * 12.;

        particles.push(Particle {
            pos: vec2(x, y),
            vel: vec2(0., 0.),
        })
    }

    particles
}

fn b2f(b: bool) -> f32 {
    if b { 1. } else { 0. }
}

fn compute_rest_density(grid: &mut MacGrid) {
    // Average density of water cells
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

fn solve_incompressibility(grid: &mut MacGrid) {
    let nx = grid.nx;
    let ny = grid.ny;

    let is_wall = |i, j| grid.cell_type[i + j * nx] == CellType::Wall;
    let iu = |i, j| i + j * (nx + 1);
    let iv = |i, j| i + j * nx;

    for _ in 0..1 {
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
                    * (grid.u[grid.iu(i + 1, j)] - grid.u[grid.iu(i, j)]
                        + grid.v[grid.iv(i, j + 1)]
                        - grid.v[grid.iv(i, j)]);

                d -= DENSITY_STIFFNESS * (grid.density[i + j * nx] - grid.rest_density);

                grid.u[iu(i, j)] += d * s_l / s;
                grid.u[iu(i + 1, j)] -= d * s_r / s;
                grid.v[iv(i, j)] += d * s_t / s;
                grid.v[iv(i, j + 1)] -= d * s_b / s;
            }
        }
    }
}

fn simulate_particles(particles: &mut [Particle], w: f32, h: f32) {
    let nx = (w / CELL_SIZE).ceil() as usize;
    let ny = (h / CELL_SIZE).ceil() as usize;

    for p in particles {
        p.vel += GRAVITY * DT;
        p.pos += p.vel * DT;

        // Wall collision detection: ensure particle center is in a non-wall cell
        // Wall cells are at i=0, i=nx-1, j=0, j=ny-1
        let min_water_x = CELL_SIZE + RADIUS;
        let max_water_x = (nx - 1) as f32 * CELL_SIZE - RADIUS;
        let min_water_y = CELL_SIZE + RADIUS;
        let max_water_y = (ny - 2) as f32 * CELL_SIZE - RADIUS;

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

#[macroquad::main("FLIP")]
async fn main() {
    let mut w = screen_width();
    let mut h = screen_height();

    let mut particles = spawn_particles(10000);
    let mut grid = MacGrid::new(w, h);

    loop {
        clear_background(BLACK);

        if w != screen_width() || h != screen_height() {
            w = screen_width();
            h = screen_height();
            grid = MacGrid::new(w, h);
        } else {
            grid.clear();
        }

        simulate_particles(&mut particles, w, h);

        particles_to_grid(&particles, &mut grid);

        compute_particle_density(&particles, &mut grid);

        compute_rest_density(&mut grid);

        // save pre-pressure grid velocities for FLIP update
        let u_prev = grid.u.clone();
        let v_prev = grid.v.clone();

        solve_incompressibility(&mut grid);

        grid_to_particles(&mut particles, &grid, &u_prev, &v_prev);

        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let cell_type = &grid.cell_type[i + j * grid.nx];
                let color = match cell_type {
                    CellType::Air => Color::new(0.1, 0.1, 0.1, 1.0),
                    CellType::Water => Color::new(0.0, 0.5, 1.0, 1.0),
                    CellType::Wall => Color::new(0.5, 0.5, 0.5, 1.0),
                };
                let x = i as f32 * CELL_SIZE;
                let y = j as f32 * CELL_SIZE;
                draw_rectangle(x, y, CELL_SIZE, CELL_SIZE, color);
            }
        }

        for p in &mut particles {
            draw_circle(p.pos.x, p.pos.y, RADIUS, SKYBLUE);
        }
        next_frame().await
    }
}
