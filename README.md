<div align="center">
<div class="header">
  <h1>PolarBearel</h1>
 [Read our paper](https://doi.org/...) | [Search the PolarBearal database](https://isitabarrel.ku.edu/search) | [Install from source](https://github.com/SluskyLab/PolarBearel3) <br> <br>
</div>
</div>

----

# PolarBearel3

A highly accurate strand assignment and barrel analysis program.

This repository contains the algorithm, source code and accompanying scripts for PolarBearel3.

If you find our work useful, please cite our accompanying paper: https://doi.org/...

## Installation

PolarBearel requires a C# runtime to build the executable.

- Windows/MacOS uses `.NET 6.0` or higher.
- Build-compatible with `Mono` latest on Linux.

To install PolarBearel from source:

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