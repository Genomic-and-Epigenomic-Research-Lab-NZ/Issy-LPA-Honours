# envs

Conda environments used in this project (Section 2.4, Tables 1-3).

| File | Used for |
|---|---|
| `custom_env.yml` | Main analysis environment (Table 1) |
| `vacmap_env.yml` | VACmap realignment (Table 2) |
| `bandage_env.yml` | Bandage assembly graph images (Table 3) |
| `samtools_remote.yml` | samtools with remote file access, for downloading the validation data |

Recreate an environment with:

    conda env create -f custom_env.yml

VACmap itself was installed from its GitHub source into `vacmap_env` and is not
listed in `vacmap_env.yml`. Dorado and KILDA were installed separately (see the
main README).
