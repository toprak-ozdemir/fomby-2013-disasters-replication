# Original replication package

The files in this folder come from the replication package of

> Fomby, T., Ikeda, Y., & Loayza, N. V. (2013). The Growth Aftermath of Natural Disasters.
> *Journal of Applied Econometrics*, 28(3), 412–434. https://doi.org/10.1002/jae.1273

The package was downloaded from the **JAE Data Archive**:
http://qed.econ.queensu.ca/jae/2013-v28.3/fomby-ikeda-loayza/

All rights stay with the original authors. The files are included only so the Python replication
in this repository runs without extra steps, and so readers can compare it with the original code.
The MIT License covering the rest of this repository does **not** apply here.
For the authoritative version, please download the package from the archive.

## Contents

| Path              | Source in the archive      | Notes |
|-------------------|----------------------------|-------|
| `readme.fil.txt`  | `readme.fil.txt`           | Unchanged. Describes the data sources and variables. |
| `matlab/*.m`      | `matlab-files.zip`         | Unchanged (byte-identical to the archive). |
| `stata/*.do`      | `stata-text-files.zip`     | Unchanged (byte-identical to the archive). |
| `data/*.csv`      | `spreadsheets-xls.zip`     | Converted from `.xls` to `.csv`, see below. |

### About the CSV files

The six estimation datasets (`dev_grow`, `dev_agr`, `dev_non`, `rich_grow`, `rich_agr`,
`rich_non`) are distributed in the archive as `.xls` (and as ASCII `.txt`). They were converted
to headerless `.csv` because the Python code reads CSV. The values were checked against the
authors' ASCII versions in `spreadsheets-text.zip` and are numerically identical: same shape,
maximum absolute difference 0.

### Not included

To keep the repository small, the following archive files are **not** copied here. Get them from
the archive if you need them:

- `stata-data-file.zip` (`fil-data.dta`) and `fil-data.txt`: the raw country-year panel used by
  `stata/data_prep.do`
- the `.xls` and `.txt` versions of the six estimation datasets

The paper and its online supplement are copyrighted by John Wiley & Sons and are also not
included. They're available through the DOI above.
