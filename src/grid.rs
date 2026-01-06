use crate::constants::CELL_SIZE;

#[derive(Debug, Clone, PartialEq)]
pub enum CellType {
    Air,
    Water,
    Wall,
}

pub struct MacGrid {
    pub nx: usize,
    pub ny: usize,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub wu: Vec<f32>,
    pub wv: Vec<f32>,
    pub cell_type: Vec<CellType>,
    pub density: Vec<f32>,
    pub rest_density: f32,
}

impl MacGrid {
    pub fn new(w: f32, h: f32) -> Self {
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

    pub fn clear(&mut self) {
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

#[inline]
pub fn b2f(b: bool) -> f32 {
    if b { 1. } else { 0. }
}
