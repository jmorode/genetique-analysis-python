# Genetics Package

This folder provides tools for analyzing genetic data, computing population genetics statistics, and performing various genetic analyses including pairwise comparisons, ML-Relate analysis, and simulation studies.

## Features

- **Population genetics analysis**: Compute allele frequencies, heterozygosity, and population genetic statistics
- **Pairwise distance calculations**: Analyze genetic distances between individuals and populations
- **ML-Relate integration**: Generate input files and analyze ML-Relate output for relationship inference
- **Family simulation**: Generate simulated families for validation and testing
- **Error analysis**: Test the impact of genotyping errors on population genetic analyses
- **Recapture analysis**: Identify potential recaptured individuals based on genetic similarity
- **Statistical testing**: Perform Kolmogorov-Smirnov tests and other statistical analyses
- **Visualization**: Generate plots and visualizations of genetic data and analysis results


## Project Structure

The package is organized into several modules:

```
src/genetics/
├── core/           # Core functionality and configuration
├── io/             # Input/output utilities
├── analysis/       # Analysis modules
└── visualization/  # Plotting utilities (future)
```

### Core Modules

- **`config.py`**: Configuration management and data loading
- **`constants.py`**: Constants and reference data used throughout the package

### Analysis Modules

- **`description.py`**: Population genetics descriptive statistics
- **`pairwise_differences.py`**: Pairwise genetic distance calculations
- **`generation.py`**: Family simulation and generation tools
- **`ml_relate.py`**: ML-Relate integration and analysis
- **`recaptures.py`**: Recapture identification
- **`test_stats.py`**: Statistical testing functions
- **`error_generation.py`**: Error simulation and impact analysis

## Required Input Files

The genetics package requires specific input files organized in a project directory structure:

```
projects/
└── {project_name}/
    └── inputs/
        ├── genotypes.csv        # Main genotype data (required)
        ├── selection_{name}.csv # Population/year selection (required)
        ├── conversion.csv       # Allele conversion table (optional)
        └── config.json         # Analysis configuration (optional)
```

### Input File Formats

#### 1. genotypes.csv (Required)
Contains the genetic data with the following columns:
- `Sample`: Individual identifier
- `Population`: Population name
- `Year`: Sampling year
- `Subcategory`: Additional categorization (optional)
- Locus columns: One column per genetic locus

Example:
```csv
Sample;Population;Year;Subcategory;Locus1;Locus2;Locus3
IND001;Pop_A;2020;Adult;120;135;200
IND001;Pop_A;2020;Adult;124;139;204
IND002;Pop_B;2021;Juvenile;118;135;196
IND002;Pop_B;2021;Juvenile;122;131;200
```

*Note: Genotypes are stored in two-line format (one line per allele)*

#### 2. selection_{name}.csv (Required)
Specifies which populations and years to include in the analysis:
```csv
Population;Year;Subcategory
Pop_A;2020;Adult
Pop_B;2021;Juvenile
Pop_C;2020;
```

#### 3. conversion.csv (Optional)
Maps raw allele codes to standardized numeric values:
```csv
Loci;allele_raw;allele_corrected
Locus1;A;100
Locus1;B;104
Locus2;X;200
Locus2;Y;210
```

#### 4. config.json (Optional)
Contains analysis-specific configurations:
```json
{
  "analysis_param1": true,
  "analysis_param2": 0.05,
  "output_options": ["plots", "tables"]
}
```

An example of a config.json file is provided in `config.example.json`.

## Aggregation Types

The package supports different aggregation strategies:

- **`"all"`**: Combine all populations and years together
- **`"pops"`**: Group by population (merge years within populations)
- **`"pop_years"`**: Separate analysis for each population-year combination
- **`"subcat"`**: Include subcategory information in grouping

## Command Line Usage

The package provides two command-line entry points:

### Main Analysis Pipeline

```bash
genetics-main
```

### Error Testing Pipeline
```bash
genetics-error-tests
```

## Analysis Workflow

The typical analysis workflow includes:

1. **Data Loading and Preprocessing**
   - Load genotype data
   - Apply allele conversions (if specified)
   - Convert two-line genotypes to single-line format

2. **Descriptive Analysis**
   - Calculate allele frequencies
   - Compute heterozygosity statistics
   - Generate population genetic summaries

3. **Pairwise Analysis**
   - Calculate genetic distances between all individual pairs
   - Generate pairwise distance matrices
   - Create distance plots and visualizations

4. **Recapture Analysis**
   - Identify potential recaptured individuals
   - Generate recapture probability matrices

5. **ML-Relate Analysis**
   - Generate ML-Relate input files
   - Process ML-Relate output
   - Validate relationship predictions with simulated families

6. **Statistical Testing**
   - Perform hypothesis tests on genetic distances
   - Compare observed vs. simulated data
   - Generate statistical reports


    
## Output Structure

Analysis results are saved in the project output directory:

```
projects/
└── {project_name}/
    └── outputs/
        ├── raw_data/                # Processed input data
        ├── heterozygosity/          # Population genetic statistics
        ├── pairwise_differences/    # Distance matrices and plots
        ├── recaptures/             # Recapture analysis results
        ├── ml_relate/              # ML-Relate files and results
        ├── test_stats/             # Statistical test results
        └── error_recapture_tests/  # Error analysis results
```

## Dependencies

- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing
- **matplotlib**: Plotting and visualization
- **scipy**: Scientific computing and statistics
- **joblib**: Parallel computing
- **tqdm**: Progress bars


## License

This project is licensed under the CeCILL free software license agreement - see the LICENSE file for details.
