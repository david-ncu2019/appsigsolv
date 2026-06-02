# appsigsolv Architecture Plan
Refactor and merge `gps_decompose.py`, `reconstruct_signal.py`, and the `omt_ncu` module into a unified, modular CLI Python package named `appsigsolv` (Applied Geology's signal solver).

## Architecture

```
appsigsolv/
├── __init__.py
├── __main__.py               (Entry point: enables `python -m appsigsolv`)
├── cli/
│   ├── __init__.py
│   ├── parser.py             (Argparse setup and subcommand routing)
│   ├── cmd_decompose.py      (Workflow logic for the 'decompose' command)
│   └── cmd_reconstruct.py    (Workflow logic for the 'reconstruct' command)
├── core/
│   ├── __init__.py
│   ├── modeling.py           (Design matrix generation, time functions, model fitting)
│   └── dia.py                (OMT calculations, Lomb-Scargle periodograms, jump detection)
├── io/
│   ├── __init__.py
│   └── data_manager.py       (Data loading, cleaning, datetime parsing, and saving outputs)
└── utils/
    ├── __init__.py
    └── visualization.py      (Matplotlib plotting and report generation)
```
