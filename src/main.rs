use macroquad::{miniquad::window::screen_size, prelude::*};

mod constants;
mod grid;
mod particle;
mod renderer;
mod simulation;

use grid::MacGrid;
use particle::{
    advect_particles, handle_particle_collisions, push_particles_apart, spawn_particles,
};
use renderer::render;
use simulation::{
    apply_gravity, compute_particle_density, grid_to_particles, particles_to_grid,
    solve_incompressibility,
};

#[cfg(target_arch = "wasm32")]
unsafe extern "C" {
    fn get_accel_x() -> f64;
    fn get_accel_y() -> f64;
    fn get_accel_z() -> f64;
}

fn get_acceleration() -> Vec2 {
    let acc = {
        #[cfg(target_arch = "wasm32")]
        unsafe {
            vec2(get_accel_x() as f32, get_accel_y() as f32)
        }

        #[cfg(not(target_arch = "wasm32"))]
        {
            vec2(0., 10.)
        }
    };

    acc * 5000.
}

#[macroquad::main("FLIP")]
async fn main() {
    let (mut w, mut h) = screen_size();

    let mut particles = spawn_particles(10000);
    let mut grid = MacGrid::new(w, h);

    loop {
        let new_size = screen_size();
        if (w, h) != new_size {
            (w, h) = new_size;
            grid = MacGrid::new(w, h);
        } else {
            grid.clear();
        }
        grid.clear();

        advect_particles(&mut particles);
        push_particles_apart(&mut particles);
        handle_particle_collisions(&mut particles, &grid);

        particles_to_grid(&particles, &mut grid);
        apply_gravity(&mut grid);

        compute_particle_density(&particles, &mut grid);

        let u_prev = grid.u.clone();
        let v_prev = grid.v.clone();

        solve_incompressibility(&mut grid);

        grid_to_particles(&mut particles, &grid, &u_prev, &v_prev);

        render(&particles);

        next_frame().await;
    }
}
