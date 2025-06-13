# The conditional saddlepoint approximation for fast and accurate large-scale hypothesis testing

This repository reproduces the results reported in the following paper:

Z. Niu, J. Ray Choudhury, E. Katsevich. "The conditional saddlepoint approximation for fast and accurate large-scale hypothesis testing."
([arXiv v2](https://arxiv.org/abs/2407.08911v2))

There are two options for reproducing results: [plotting from downloaded results](#option-1-quick-figures-from-downloaded-results) or [complete reanalysis](#option-2-complete-rerun-all-analyses).

## Option 1 (quick): Figures from Downloaded Results

**Requirements:** Laptop or desktop with R installed and a Unix-like OS (e.g. MacOS or Windows WSL2). Package dependencies are automatically installed by the scripts via `renv`.

### Steps

1. **Clone repository**
   ```bash
   git clone git@github.com:Katsevich-Lab/spacrt-manuscript.git
   cd spacrt-manuscript
   ```

2. **Download results** from [Dropbox](https://www.dropbox.com/scl/fo/02r52y06g7p6h378ojsuk/AMPvwTqGBHWxnaZqqZfel9I?rlkey=9ta49vbw966rfljng3e2p261x&st=s2xjfxf4&dl=0).

3. **Configure paths**
   ```bash
   echo 'LOCAL_SPACRT_DATA_DIR="/path/to/downloaded/results"' > ~/.research_config
   ```

4. **Generate figures and tables**
   ```bash
   ./run_all_simulation.sh plot-only
   ./run_all_realdata.sh plot-only
   ```

Figures and tables will be saved to `manuscript/figures-and-tables/`.

## Option 2 (Complete): Rerun All Analyses

**Requirements:** HPC cluster with SGE job scheduler with R installed. Package dependencies are automatically installed by the scripts via `renv`.


### Steps

1. **Clone repository**
   ```bash
   git clone git@github.com:Katsevich-Lab/spacrt-manuscript.git
   cd spacrt-manuscript
   ```

2. **[Install](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [configure](https://www.nextflow.io/docs/latest/config.html) Nextflow for your cluster**

3. **Download Gasperini dataset** from [Dropbox](https://www.dropbox.com/scl/fo/sc7c5ezm45piaicyi08ia/AImzh7IgxoyzS68JAx8HOLA?rlkey=avts6iadwq9zdhmy8t2kyjb0v&st=cw8cdsj6&dl=0). These data were obtained from the [published data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861) using [this pipeline](https://github.com/Katsevich-Lab/import-gasperini-2019-v2).

4. **Configure paths**
   ```bash
   echo 'LOCAL_SPACRT_DATA_DIR="/path/for/output/results"' > ~/.research_config
   echo 'LOCAL_GASPERINI_2019_V2_DATA_DIR="/path/to/gasperini/data"' >> ~/.research_config
   ```

5. **Run complete analysis**
   ```bash
   qsub run_all_simulation.sh
   qsub run_all_realdata.sh
   ```

Results will be saved to `/path/for/output/results`; figures and tables will be saved to `manuscript/figures-and-tables/`.

</details>
