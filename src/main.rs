use macroquad::{miniquad::window::screen_size, prelude::*};

mod constants;
mod grid;
mod particle;
mod renderer;
mod simulation;

use grid::MacGrid;
use particle::{
    handle_particle_collisions, integrate_particles, push_particles_apart, spawn_particles,
};
use renderer::render;
use simulation::{
    compute_particle_density, grid_to_particles, particles_to_grid, solve_incompressibility,
};

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
