# FLIP Fluid Simulation

A fluid simulation using the FLIP (Fluid-Implicit Particle) method implemented in Rust with macroquad.

## Local Development

### Prerequisites

- Rust with wasm32-unknown-unknown target: `rustup target add wasm32-unknown-unknown`
- Python 3 (for local server)

### Building for Web

```bash
# Build WASM for deployment
./build.sh

# Build and serve locally for testing
./serve.sh
```

## Deployment

This project is automatically deployed to GitHub Pages when pushed to the main branch.

### Deployment Structure

- `dist/` - Built assets (tracked for GitHub Pages)
- `target/` - Rust build artifacts (ignored)
- `.github/workflows/deploy.yml` - CI/CD pipeline

The deployment:

1. Builds the Rust project for `wasm32-unknown-unknown` target
2. Copies the generated `.wasm` file to `dist/`
3. Updates HTML to reference local paths
4. Deploys to GitHub Pages using official GitHub Actions

## Technical Details

- **Grid size**: Adaptive based on window dimensions (10px cells)
- **Particles**: 10,000 water particles
- **Simulation steps**:
  - Particle integration with gravity
  - Particle separation (3 iterations)
  - Grid transfer (particles → grid)
  - Density computation
  - Incompressibility solver (5 iterations)
  - Velocity correction (grid → particles)
- **Collision handling**: Boundary box walls

## License

MIT License

