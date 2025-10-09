# Genetique Analysis

An analysis package for microsatellites genetic data processing, providing tools for frequency calculations, pairwise analysis, ML-related analysis, and statistical testing.

## Features

- **Analysis pipeline for**:
  - Frequency and heterozygosity calculations
  - Pairwise difference analysis and plotting
  - Recapture analysis
  - ML-related analysis and reliability testing
  - Statistical testing and relationship analysis

## Installation

### Using uv (recommended)

```bash
# Install from local source
uv pip install -e .
```

### Using pip

```bash
# Install from local source
pip install -e .
```

### Using the provided install script

```bash
# The script automatically detects and uses uv or pip
./install.sh
```

## Quick Start

1. **Create a project directory structure**:
   ```
   projects/
   └── your_project/
       ├── inputs/
       │   ├── config.json
       │   ├── genotypes.csv
       │   ├── selection_{selection_name}.csv
       │   └── conversion.csv (optional)
       └── outputs/
   ```

2. **Create a configuration file** (`projects/your_project/inputs/config.json`):
   ```json
   {
     "project_name": "your_project",
     "selection_name": "everything",
     "aggregation_type": "pops",
     "population_order": [],
     "analysis_steps": {
       "step_1_frequency_heterozygosity": true,
       "step_2_pairwise_differences": true,
       "step_3_plot_pairwise_distances": true,
       "step_4_recaptures": true,
       "step_5_ml_relate_input": false,
       "step_6_ml_relate_output": false,
       "step_7_relationship_analysis": false
     },
     "simulation_parameters": {
       "nb_families": 1000,
       "nb_individuals": 100
     }
   }
   ```

3. **Run the analysis**:
   ```bash
   # Using the command-line interface
   genetique-analysis your_project
   
   
   # Or using Python directly
   python -m genetique_analysis.cli.main your_project
   ```

## Command Line Usage

### Main Analysis Pipeline

```bash
# Run analysis for a specific project
genetique-analysis my_project_name

# Run analysis for the default project
genetique-analysis

```

### Error Impact Testing

```bash
# Run error impact tests with default parameters
genetique-error-tests

# Run with custom parameters
genetique-error-tests project_name selection_name population_name nb_individuals nb_iterations

# Example with specific parameters
genetique-error-tests analysis_nov_2024_more_loci pachon_sabinos CF_Pachon 100 1000
```

### Alternative Python Module Usage

```bash
# Run main analysis
python -m genetique_analysis.cli.main your_project_name

# Run error tests
python -m genetique_analysis.cli.error_tests project_name selection_name pop_name nb_inds nb_iterations
```

## Configuration System

The main functions of the package are configuration-driven. All analysis parameters and steps are controlled through a JSON configuration file.

### Configuration File Location

Place your configuration file at: `projects/{project_name}/inputs/config.json`

### Configuration File Structure

```json
{
  "project_name": "your_project_name",
  "selection_name": "your_selection_name",
  "aggregation_type": "pops",
  "population_order": ["Pop1", "Pop2", "Pop3"],
  "analysis_steps": {
    "step_1_frequency_heterozygosity": true,
    "step_2_pairwise_differences": true,
    "step_3_plot_pairwise_distances": true,
    "step_4_recaptures": true,
    "step_5_ml_relate_input": true,
    "step_6_ml_relate_output": true,
    "step_7_relationship_analysis": true
  },
  "simulation_parameters": {
    "nb_families": 1000,
    "nb_individuals": 100
  },
  "plotting_parameters": {
    "figure_size": [8, 6],
    "large_figure_size": [30, 10],
    "dpi": 300
  },
  "file_parameters": {
    "csv_separator": ";",
    "file_encoding": "utf-8"
  }
}
```

### Configuration Parameters

#### Basic Parameters
- **project_name**: Name of your project directory
- **selection_name**: Name of the selection file (without .csv extension)
- **aggregation_type**: Analysis aggregation type
  - `"all"`: All populations together
  - `"pops"`: Each population separately (merge years)
  - `"pop_years"`: Each population/year combination separately
  - `"subcat"`: Account for subcategories
- **population_order**: List of population names in desired order for plotting (empty list uses default order)

#### Analysis Steps Control
You can enable/disable individual analysis steps:

- **step_1_frequency_heterozygosity**: Frequency and heterozygosity calculations
- **step_2_pairwise_differences**: Pairwise difference calculations
- **step_3_plot_pairwise_distances**: Pairwise distance plotting
- **step_4_recaptures**: Recapture analysis
- **step_5_ml_relate_input**: ML-related input file generation
- **step_6_ml_relate_output**: ML-related output analysis
- **step_7_relationship_analysis**: Relationship analysis and statistical testing

#### Simulation Parameters
- **nb_families**: Number of families for simulation (default: 1000)
- **nb_individuals**: Number of individuals for simulation (default: 100)

#### Plotting Parameters
- **figure_size**: Default figure size [width, height] (default: [8, 6])
- **large_figure_size**: Large figure size [width, height] (default: [30, 10])
- **dpi**: Figure resolution (default: 300)

#### File Parameters
- **csv_separator**: CSV file separator (default: ";")
- **file_encoding**: File encoding (default: "utf-8")

### Example Configurations

#### 1. Full Analysis (All Steps Enabled)
```json
{
  "project_name": "my_project",
  "selection_name": "everything",
  "aggregation_type": "pops",
  "population_order": [],
  "analysis_steps": {
    "step_1_frequency_heterozygosity": true,
    "step_2_pairwise_differences": true,
    "step_3_plot_pairwise_distances": true,
    "step_4_recaptures": true,
    "step_5_ml_relate_input": true,
    "step_6_ml_relate_output": true,
    "step_7_relationship_analysis": true
  }
}
```

#### 2. Skip ML-Related Analysis (Steps 5, 6, 7 Disabled)
```json
{
  "project_name": "my_project",
  "selection_name": "everything",
  "aggregation_type": "pops",
  "population_order": [],
  "analysis_steps": {
    "step_1_frequency_heterozygosity": true,
    "step_2_pairwise_differences": true,
    "step_3_plot_pairwise_distances": true,
    "step_4_recaptures": true,
    "step_5_ml_relate_input": false,
    "step_6_ml_relate_output": false,
    "step_7_relationship_analysis": false
  }
}
```

### Configuration Requirements

- **Configuration file is required**: The system requires a `config.json` file to run
- **Missing parameters**: If a configuration file exists but is missing some parameters, default values are used
- **Error handling**: Clear error messages are provided if the configuration file is missing or invalid

## Python API

You can also use the package programmatically:

### Main Analysis Pipeline

```python
from genetique_analysis import FileConfiguration
from genetique_analysis.analysis import (
    compute_frequency_and_heterozygosity,
    save_proba_same_ind,
    recover_pairwise_difference,
    get_recaptures
)

# Load configuration
config = FileConfiguration("my_project", "pops", "everything")
config.initialize_global()

# Run specific analysis steps
save_proba_same_ind(config)
compute_frequency_and_heterozygosity(config, [])
recover_pairwise_difference(config)
get_recaptures(config)
```

### Error Impact Testing

```python
from genetique_analysis.cli import main_error_recaptures

# Run error impact testing
main_error_recaptures(
    project_name="my_project",
    selection_name="everything", 
    pop_name="Population1",
    nb_inds_simu=100,
    nb_iterations_simu=1000
)
```

## Project Structure

```
genetique_analysis/
├── __init__.py              # Main package with lazy imports
├── core/                    # Core functionality
│   ├── __init__.py         # Core module with lazy imports
│   ├── config.py           # Configuration management
│   └── constants.py        # Constants and defaults
├── utils/                   # Utility functions
│   ├── __init__.py         # Utils module exports
│   └── utils.py            # Helper functions
├── analysis/                # Analysis modules
│   ├── __init__.py         # Analysis module with lazy imports
│   ├── description.py      # Frequency and heterozygosity
│   ├── pairwise_differences.py  # Pairwise analysis
│   ├── ml_relate.py        # ML-related analysis
│   ├── recaptures.py       # Recapture analysis
│   ├── generation.py       # Simulation and generation
│   ├── error_generation.py # Error simulation
│   └── test_stats.py       # Statistical testing
└── cli/                     # Command-line interface
    ├── __init__.py         # CLI module exports
    └── main.py             # CLI entry point
```

## Requirements

- Python 3.9+
- NumPy
- Pandas
- Matplotlib
- SciPy
- Loguru
- Joblib
- tqdm


## License

This project is licensed under CeCILL free software license agreement - see the LICENSE file for details.
