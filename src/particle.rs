use macroquad::prelude::*;
use std::collections::HashMap;

use crate::constants::{CELL_SIZE, DT, RADIUS, SEPARATION_STEPS};
use crate::grid::MacGrid;

pub struct Particle {
    pub pos: Vec2,
    pub vel: Vec2,
}

pub fn integrate_particles(particles: &mut [Particle]) {
    for p in particles {
        p.vel += crate::get_acceleration() * DT;
        p.pos += p.vel * DT;
    }
}

pub fn push_particles_apart(particles: &mut [Particle]) {
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

pub fn handle_particle_collisions(particles: &mut [Particle], grid: &MacGrid) {
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

pub fn spawn_particles(n: usize) -> Vec<Particle> {
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
