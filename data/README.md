# data/

`Rfam.cm` and `Rfam.clanin` sourced 03/08/2026 from [rfam.org](https://rfam.org) ([https://ftp.ebi.ac.uk/pub/databases/Rfam/](https://ftp.ebi.ac.uk/pub/databases/Rfam/)), Rfam v15.1.

```
wget -O Rfam.cm https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.cm.gz && gunzip Rfam.cm.gz
wget -O Rfam.clanin https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.clanin
```

`Rfam.cm` can be passed directly to either tool's `-cm-file` flag; every named model inside is expanded automatically via `cmstat`/`cmfetch` (see the top-level README).
