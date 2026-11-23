# CAD sample tasks

Nine engineering tasks for coding agents, in Harbor format. Each task gives the agent
real CAD or engineering files and asks for a working program or model. A hidden test
suite scores the result 0 or 1.

| Task | What the agent builds |
|---|---|
| `gbroques-openafpm-cad-core-e0767e2-b173cc85` | Parametric OpenAFPM Furling Clearance and Stop Design |
| `pachterlab-poseidon-pump-carriage-v2-7b8fb975` | Create Print-Fit Poseidon Hardware with 0.40 mm Bore Clearance |
| `pachterlab-poseidon-pump-carriage-v2-8aa8f3a5` | Implement STL Bore Access Survey Tool |
| `stuttgart-timber-beam-columns-f1589626` | Generate and Rate a Glulam Beam-Column Layup Release |
| `tum-gni-anonymized-course-build-d0a802ff-0585e632` | IFC Truss Load-Rating and Operating Chart Generator |
| `vendor-cad-assemblies-4ab00b55-d4fc6754` | Unfold Folded Sheet-Metal STLs into Laser Blanks with Brake Marks |
| `vendor-cad-assemblies-758f56a9-91411172` | Generate Mold-Pull Trapped-Steel Solids from STL Shells |
| `zenodo-krijnen-balcony-ifc-models-a8-51eb5f01` | IFC Cast-in-Place Concrete Takeoff Program |
| `zipcpu-wbscope-85907874` | Implement Windowed DDR Capture Reassignment in SystemVerilog |

Measured results for three models are in each task's `task.toml`.

## Run

```sh
harbor run --repo Scalable-Worlds/cad-samples -d cad-samples -a <agent> -m <model>
```

From a clone, one task at a time:

```sh
harbor run -p tasks/<task> -a oracle    # the reference solution scores 1
```

Each task pins a prebuilt image in a private registry. For pull access or any
question, write to contact@scalableworlds.ai.

Shared privately. Please do not redistribute.
