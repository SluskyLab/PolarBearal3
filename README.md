<div align="center">
<div class="header">
  <h1>PolarBearel 3.0</h1>
  <p>A highly accurate strand assignment and barrel analysis program.<p>
</div>

  [Read our paper on cell.com](https://www.cell.com/biophysj/fulltext/S0006-3495(26)00188-8) | [Install from source](https://github.com/SluskyLab/PolarBearel3) | [Search the PolarBearal database](https://isitabarrel.ku.edu/search) <br> <br>
</div>

Introduction
----

This repository contains the algorithm, source code, and accompanying scripts for the latest version of PolarBearel3.

If you find our work useful, please consider citing our accompanying paper now published in the Biophysical Journal (Cell Press): https://doi.org/10.1016/j.bpj.2026.03.016

## Installation

PolarBearel requires a C# runtime to build the executable.

- Windows/MacOS uses `.NET 6.0` or higher.
- Build-compatible with `Mono` latest on Linux.

### To install PolarBearel from source:

1. Clone the repository from Github.

```sh
git clone https://github.com/SluskyLab/PolarBearel3
```

2. Build the PolarBearel executable.
```sh
cd PolarBearel3 && dotnet build
```

## Reproducing PolarBearel strand assignments from AlphaFoldDB

> Note: PolarBearelDB strand number assignments are already available for [public download](https://isitabarrel.ku.edu/download) as part of our [IsItABarrel web application](https://isitabarrel.ku.edu/).

Download the AlphaFold structures necessary for barrel and strand counting analysis.

```sh
cd data/
chmod +x download_afdb_structures.sh
./download_afdb_structures.sh
```

After the structures have been downloaded, follow the menu to set your current dataset to the AlphaFold structures. The process has been automated through the example file `commands.txt`.

```sh
dotnet run < commands.txt
```

Given the large size of the dataset, we recommend running this as a detached or background process to avoid incomplete results.

### Accessing predictions and results

Each PDB structure receives its own strand assignment and annotated PyMOL session. Batch results can be found under `Output/...`.

## License

PolarBearel is provided under the GPL 3.0 Public License. For more information, refer to our [LICENSE](./LICENSE.md) file.
