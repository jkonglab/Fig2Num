# Introduction

We develop a tool

Prediction of outbreaks and spatial spread of infectious diseases highly depends on accurately collecting and integrating large-scale epidemic data from multiple resources. However, such data are often demonstrated in bar chart figures published in a format where original data can be readily used. As a result, the capability to extract such data from bar charts in an automated fashion becomes the key to enabling ensuing prediction analysis. We develop tool for efficient and accurate data extraction from bar chart figures. The developed conversion method leverages the power of image analysis and only requires minor manual intervention, enabling analysis in a batch mode for large-scale data conversion.


# Usage

(1)
```bash
# extract data from a bar chart in style 1:
figure2num_clean(data_path, filename, result_path)
```
e.g.
```bash
figure2num_clean('../data/', 'figure1_20190808.png', '../result/')
```

(2)
```bash
# extract data from a bar chart in style 2:
figure2num_DRC(data_path, filename, dates, result_path)
```
e.g.
```bash
figure2num_DRC('../data/', 'SITREP_EVD_DRC_20191126-eng.png', [2018 4 30 2019 11 18], '../result/')
```

(3)
```bash
# extract data from a bar chart in style 3 (i.e. in grouped style):
figure2num_group(data_path, filename, dates, result_path)
```
e.g.
```bash
figure2num_group('../data/', 'nejmsr_fig4.tiff', [2018 8 1 2019 5 1], '../result/')
```

# Programming Language
Matlab
