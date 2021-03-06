# position-based-fluids

[![macOS](https://github.com/yuki-koyama/position-based-fluids/actions/workflows/macos.yml/badge.svg)](https://github.com/yuki-koyama/position-based-fluids/actions/workflows/macos.yml)

An implementation of Position Based Fluids [Macklin+, SIGGRAPH 2013] in C++.

![](./docs/sample-1.gif)
![](./docs/sample-2.gif)

## Dependencies

- [alembic](https://github.com/alembic/alembic) (included as a git submodule)
- [Eigen](https://eigen.tuxfamily.org/)
- [parallel-util](https://github.com/yuki-koyama/parallel-util) (included as a git submodule)
- [timer](https://github.com/yuki-koyama/timer) (included as a git submodule)

## Prerequisites

macOS:
```sh
brew install cmake eigen imath
```

Other environments (e.g., Ubuntu, Windows) are not tested.

## Implementation Details

Supported:
- Incompressibility constraint
- Artificial surface tension
- XSPH viscosity
- Hash-grid neighbor search
- Parallelization based on multi-threading
- Alembic export

Not supported:
- Vorticity confinement
- GPU computing
- Surface reconstruction
- Rendering

## License

MIT License
