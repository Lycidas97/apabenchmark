# APA Benchmark

This repository contains the code used in the study "Benchmarking alternative polyadenylation detection in single-cell and spatial transcriptomes".

## Overview

Alternative polyadenylation (APA) is a widespread post-transcriptional regulatory mechanism that generates mRNA isoforms with diverse 3' UTRs and terminal exons. This study systematically compares seven 3'-tag-based sequencing protocols and six computational tools for their ability to detect APA events in single-cell and spatial transcriptomes.

## Contents

- `data/`: Contains raw data, intermediata data and result data of the study.
- `scripts/`: Contains the scripts used for data processing, analysis, and evaluation of APA detection tools. Each section within the `scripts/` directory has its own README file with detailed instructions.

## Usage

1. Clone the repository:
   ```
   git clone https://github.com/Lycidas97/apabenchmark.git
   ```

2. Install the required dependencies as specified in the documentation.

3. Navigate to the `scripts/` directory:
   ```
   cd apabenchmark/scripts/
   ```

4. Follow the instructions provided in the README file of each section within the `scripts/` directory to run the corresponding scripts. Each section focuses on a specific aspect of the analysis, such as data preprocessing, APA detection, or evaluation.

5. Run the scripts in the order specified by the README files to reproduce the analysis and generate the results.

6. The generated output files and figures will be saved in the `data/int_data/` directory.
7. 
## Contact

For questions or feedback, please contact:

- Sida Li (sida.lycidas@foxmail.com)

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

